#include "gdal_priv.h"
#include "imad.h"
#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;
using std::cout; //debugging purposes
using std::endl;

//typedef Map<Matrix<double,Dynamic,Dynamic,RowMajor> > MapRMMatrixXd;
//defined in imad.h

void imad(std::string filename1="", std::string filename2="",
          std::string output_file="", std::string format="",
          int* bands1_arg=NULL, int* bands2_arg=NULL, int input_bands=-1,
          int* win_sz_arg=NULL, int* data1_offsets=NULL,int* data2_offsets=NULL,
          double pen_inp=0.0, int maxiterations = 1, double err_tol = 0.0001
          ){


  /* We start off with a large block of largely uninteresting error-checking
   * code. You can mostly skip all this junk, it just checks inputs for sanity*/

  GDALAllRegister(); //Must be called at start, see GDAL API for details

  /* Copy over pointer args. If any inputs are NULL, we will generate a struct
   * from user input. Later, we can free these structs by seeing if the original
   * argument was NULL. */
  int* bandlist_1 = bands1_arg;
  int* bandlist_2 = bands2_arg;
  int* win_size   = win_sz_arg;
  int* offsets_1  = data1_offsets;
  int* offsets_2  = data2_offsets;

  //openFile() will ask us for filenames if the names are blank
  GDALDataset* file1 = GdalFileIO::openFile(filename1);
  GDALDataset* file2 = GdalFileIO::openFile(filename2);
  GdalFileIO::getOutputFileInfo(output_file, format);

  //If any of the inputs are NULL, ask user to fill in values
  GdalFileIO::fix_missing_band_data(&bandlist_1, &bandlist_2, input_bands);
  GdalFileIO::fix_missing_dims_data(&win_size, &offsets_1, &offsets_2);

  //Check for any malformed inputs
  bool has_error = GdalFileIO::has_errors( file1, file2,
                        bandlist_1, bandlist_2, input_bands,
                        win_size, offsets_1, offsets_2, pen_inp);

  if(has_error){
    cout << "Error checking has determined that a parameter to imad()" <<
            "is not sane!" << endl << "Exiting now..." << endl;
    return;
  }

  int nBands = input_bands;
  int ncol = win_size[0]; //X-size is number of columns
  int nrow = win_size[1]; //Y-size is number of rows

  //Image offsets: x10 = xoffset image 1, y20 = yoffset image 2, etc.
  int x10 = offsets_1[0];
  int y10 = offsets_1[1];
  int x20 = offsets_2[0];
  int y20 = offsets_2[1];
  int maxiter = maxiterations;
  double tolerance = err_tol;
  double pen = pen_inp; //Penalty to algorithm (usually set to zero)

  //These will contain pointers to the band objects, instead of band numbers
  GDALRasterBand** bands_1  = new GDALRasterBand*[nBands];
  GDALRasterBand** bands_2  = new GDALRasterBand*[nBands];
  double* noDataValues      = new double[nBands];

  //Get handles to raster bands -- makes it easer to reference later
  for(int i = 0; i < nBands; i++){
    bands_1[i]      = (file1->GetRasterBand(bands1_arg[i]));
    bands_2[i]      = (file2->GetRasterBand(bands2_arg[i]));
    noDataValues[i] = bands_1[i]->GetNoDataValue();
  }

  //Input checking/setup is done. Prepare for calculations.

  /* We set up the "chunking" procedure for large images. In most cases, bufsize
   * will be ncol and nb_pr (number of buffers per row) will be 1. this_bufsize
   * is eq. to bufsize, but can be shrunk at the end of a row to acct. for
   * buffersize not dividing always ncol evenly (e.g. bufsize=3, ncol = 11) */

  int bufsize = -1; //size of one row in the tile buffer
  int nb_pr = imad_bigfun::find_chunksize(bufsize, ncol, nBands);
  int this_bufsize = 0;

  //Bufsize was the size for one band. We have nBands per image and two images
  double* tile = new double[2 * nBands * bufsize];

  //Initialize some values for the main calculation below
  ImageStats cpm = ImageStats(nBands * 2); //Tracks image stats (mean + cov)
  double delta = 100000; //Max change between correlations
  VectorXd rho, oldrho = VectorXd::Zero(nBands);
  MatrixXd A , B; //Eigenvector matrices (eigens in columns)
  VectorXd  sigMADs, means1, means2;
  //Declare a special vector for calculations later on
  const VectorXd one = VectorXd::Constant(nBands, 1); //col vector of all "1"

  /* Start calculations (note double conditional in for-loop: we will continue
   * until the delta is smaller than the tolerance or to maxiter iterations) */

  /* This stupid chunking algorithm makes the code utterly unreadable. The end
   * result of the loops over row, bufnum, and k is to get the mean and covar
   * of all the images, piece by piece. If iter > 0, we weight the pixels based
   * on a Chi2 statistic. All we really want is a weighted mean / covariance. */
  for(int iter = 0; iter < maxiter && delta > tolerance; iter++){
    //Start gathering image statistics

    for(int row = 0; row < nrow; row++){
      for(int bufnum = 0; bufnum < nb_pr; bufnum++){
        for(int k = 0; k < nBands; k++){
          int xstart = bufsize * bufnum;
          this_bufsize = min(bufsize, ncol - xstart);
          imad_bigfun::readToBuf((tile + k), bands_1[k],
                                 x10 + xstart, y10 + row,
                                 this_bufsize, nBands);
          imad_bigfun::readToBuf((tile + nBands + k), bands_2[k],
                                 x20 + xstart, y20 + row,
                                 this_bufsize, nBands);
        }
        if(iter > 0){ //Repeated iterations use the previous iterations as weights
          VectorXd weights(this_bufsize);
          weights = imad_bigfun::calc_weights(tile, weights, A, B, means1, means2,
                                                sigMADs, this_bufsize, nBands);
          cpm.update(tile, weights.data(), this_bufsize, 2*nBands);
        }
        else{
          cpm.update(tile, NULL, this_bufsize, 2*nBands);
        }
      }
    }
    // Finished gathering image statistics (into cpm).
    MatrixXd S = cpm.get_covar();
    VectorXd means = cpm.get_means();
    cpm.reset(); //Reset image accumulator (sets all elements to zero)

    //Set up the matrices needed as input to the eigenproblem
    MatrixXd s11 = S.block(0,0,nBands,nBands);           //Self-covariance
    MatrixXd s12 = S.block(0,nBands,nBands,nBands);      //Cross-covariance
    MatrixXd s21 = S.block(nBands,0,nBands,nBands);      //S = S.transpose()
    MatrixXd s22 = S.block(nBands,nBands,nBands,nBands); //Self-covariance
    s11 = (1 - pen) * s11  + pen * MatrixXd::Identity(nBands,nBands);
    s22 = (1 - pen) * s22  + pen * MatrixXd::Identity(nBands,nBands);
    MatrixXd b1  = s11;
    MatrixXd c1  = s12 * s22.inverse() * s21;
    MatrixXd b2  = s22;
    MatrixXd c2  = s21 * s11.inverse() * s12;

    //Solve the eigenproblem
    VectorXd mu2, mu;
    if(nBands == 1){
      mu2(0) = c1(0,0)/b1(0,0);
      A(0,0) = 1 / sqrt(b1(0,0));
      B(0,0) = 1 / sqrt(b2(0,0));
    }
    else{
      Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver1(c1,b1);
      Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver2(c2,b2);
      if (solver1.info() != Eigen::Success || solver2.info() != Eigen::Success){
        throw std::runtime_error("Eigensolver did not converge!");
      }
      mu2 = solver2.eigenvalues(); //Eigenvalues should be identical, use either
      A = solver1.eigenvectors();
      B = solver2.eigenvectors();
    }

    /* We now need to sort the eigenvals and their eigenvecs. mu2 should be
     * sorted in desc. order and eigenvecs (A and B) reordered appropriately*/

    imad_utils::reorder_eigens(mu2, A, B);

    mu = mu2.array().sqrt().matrix(); //Element-wise square root (mu2 = rho**2)

    //Calculate the variances of the eigenvectors (Var(e1), Var(e2), etc.)
    //Note that in Python, these are row vectors, while they are column in C++
    VectorXd eigvarA = (A.transpose() * A).diagonal();
    VectorXd eigvarB = (B.transpose() * B).diagonal();

    /* Calculate the penalized MAD significances. Check the README for more
     * information on the source of the variable names */

    //if pen == 0, sigma = 2*(1 - mu) and rho = mu

    VectorXd sigma = ((2*one - pen*(eigvarA + eigvarB))          //numerator
                                    /                            //---------
                (1 - pen) - 2 * mu).array().sqrt().matrix();     //denominator

    VectorXd rho = (
      mu.array() * (1 - pen)                                       //numerator
                 /                                                 //---------
((one - pen*eigvarA).array() * (one - pen*eigvarB).array()).sqrt() //denominator
                   ).matrix();

    delta = (rho - oldrho).maxCoeff(); //The max of diffs between correlations
    oldrho = rho;

    cout << "Iteration: " << iter << "  --   Delta: " << delta << endl;

    sigMADs = sigma;
    means1 = means.block(0,0,nBands,1);
    means2 = means.block(nBands,0,nBands,1);

    imad_bigfun::math_cleanup(A,B,s11,s12);

    }

    // //For debuging
    // A = MatrixXd::Zero(nBands, nBands);
    // B = MatrixXd::Zero(nBands, nBands);
    // sigMADs = VectorXd::Zero(nBands);

  //End iterations. Gear up to write final result to file.

  GDALDriver* outdriver=GetGDALDriverManager()->GetDriverByName(format.c_str());
  GDALDataset *outfile = outdriver->Create(output_file.c_str(),ncol,nrow,
                                           nBands + 1, GDT_Float64, NULL);


  GdalFileIO::writeOutputToFile(outfile, tile, A, B, x10, y10, x20, y20,
                                ncol, nrow, nb_pr, bufsize,
                                bands_1, bands_2, file1, sigMADs);

  //Cleanup in aisle 7!
  GDALClose( (GDALDatasetH) outfile);
  GDALClose( (GDALDatasetH) file1 );
  GDALClose( (GDALDatasetH) file2 );

  /* If we had to fill in any NULL arguments, the original argument pointers
   * will be null (because we made copies of them to fill in with user data).
   * Use that fact to free only the data we were forced to generate ourselves */
   if(bands1_arg == NULL) delete[] bandlist_1;
   if(bands2_arg == NULL) delete[] bandlist_2;
   if(win_sz_arg == NULL) delete[] win_size;
   if(data1_offsets == NULL) delete[] offsets_1;
   if(data2_offsets == NULL) delete[] offsets_2;
  delete[] tile;
}

int main(){ //dummy main
  std::string file1 = std::string("/home/chipbuster/POP_2014/lndsr.LE71960531999293EDC00.tif");
  std::string file2 = std::string("/home/chipbuster/POP_2014/lndsr.LE71960532000120EDC00.tif");
  std::string out = std::string("/home/chipbuster/POP_2014/C++_asd.tif");
  std::string fmt = std::string("GTiff");
  int bands1[3] = {1,2,3};
  int bands2[3] = {1,2,3};
  int win[2] = {7801, 6961};
  int offset[2] = {0,0};
  int nBands = 3;
  imad(file1,file2,out,fmt, bands1, bands2, nBands, win, offset, offset);
  return 0;
}
