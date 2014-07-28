#include "gdal_priv.h"
#include "imad.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <math.h>
#include <string.h> //Gdal libraries use const char*, need strcmp()

using namespace std;

namespace GdalFileIO{

    GDALDataset* openFile(string filename){
      /* Returns a GDALDataset handle to the file given by file1.
       * If the file cannot be opened, an exception is thrown */

      if(filename.empty()){
        cout << "Please enter a file to be opened: ";
        getline(cin, filename);
      }

      GDALDataset* dataset1 = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);

      if(dataset1 == NULL){
        throw std::invalid_argument("Error while opening file" + filename);
      }
      return dataset1;
    }

    /*************************************************************************/
    /* Checks input paramenters for all potential errors. This function should
     * never trigger in the R version because R should perform error checking */

    bool has_errors(GDALDataset* file1, GDALDataset* file2,
                   int* bands1_arg, int* bands2_arg, int nBands,
                   int* win_size, int* offsets_1, int* offsets_2, int inp_pen ){

      // Check whether input bands are within limits.


      //offsets cannot be negative
      if(offsets_1[0] < 0 || offsets_2[0] < 0 ||
         offsets_1[1] < 0 || offsets_2[1] < 0 ){
        cout << "Error: image offsets cannot be negative!" << endl;
        return true;
      }

      //Image dims can't be negative either
      if(win_size[0] < 0 || offsets_2[1] < 0){
        cout << "Error: image dimensions cannot be negative!" << endl;
        return true;
      }

      //Xoffset + Xwidth cannot exceed image Xsize
      if(offsets_1[0] + win_size[0] > file1->GetRasterXSize()){
        cout << "Error in Image 1: Dimensions too large." << endl;
        cout << "X-offset + window width exceeds image x-dimension!" << endl;
        return true;
      }

      //Perform the same check for image 2
      if(offsets_2[0] + win_size[0] > file2->GetRasterXSize()){
        cout << "Error in Image 2: Dimensions too large." << endl;
        cout << "X-offset + window width exceeds image x-dimension!" << endl;
        return true;
      }

      //Yoffset + Ywidth cannot exceed image Xsize
      if(offsets_1[1] + win_size[1] > file1->GetRasterYSize()){
        cout << "Error in Image 1: Dimensions too large." << endl;
        cout << "Y-offset + window height exceeds image y-dimension!" << endl;
        return true;
      }

      //Perform the same check for image 2
      if(offsets_2[1] + win_size[1] > file2->GetRasterYSize()){
        cout << "Error in Image 2: Dimensions too large." << endl;
        cout << "Y-offset + window height exceeds image y-dimension!" << endl;
        return true;
      }

      //Image penalty must be between 0 and 1
      if(inp_pen > 1 || inp_pen < 0){
        cout << "Error: penalization value must be between 0 and 1!" << endl;
        return true;
      }

      if(nBands < 1){
        cout << "Error: must have more than zero bands!" << endl;
        return true;
      }

      return false;
    }

    /*************************************************************************/

    int* selectBands(int nBands){
      cout << "You have specified " << nBands << " bands will be used." << endl;
      cout << "Please enter the band numbers you would like to use, separated "
           << "by a space. Press 'enter' when done." << endl;

      int* retval = new int[nBands];
      string input;
      cin >> input;
      for(int i = 0; i < nBands; i++){
        retval[i] = atof(input.c_str());
        cin >> input;
      }
      return retval;
    }

    /*************************************************************************/

    /* If explicit values are not given for bands (i.e. bands to use and number
     * of bands), this function will prompt the user to fill them */
    void fix_missing_band_data(int** bands1, int** bands2, int& nBands){
      if(nBands == -1){
        cout << "Please enter the number of bands you would like to test." << endl;
        cin >> nBands;
      }
      if (nBands < 1){
        throw std::out_of_range("Invalid number of bands.");
      }
      if(*bands1 == NULL){
        cout << "Selecting bands for image 1." << endl;
        *bands1 = selectBands(nBands);
      }
      if(*bands2 == NULL){
        cout << "Selecting bands for image 2." << endl;
        *bands2 = selectBands(nBands);
      }
    }
    /*************************************************************************/

    /* If explicit values are not given for dims (i.e. size of window and
     * image offsets), this function will prompt the user to fill them */
    void fix_missing_dims_data(int** win_size, int** offset_1,int** offset_2){
      //First, get window size
      int nrow,ncol;
      if(*win_size == NULL){
        cout << "Please enter the number of columns in the image (X-size): ";
        cin >> ncol;
        cout << "Please enter the number of rows in the image (Y-size): ";
        cin >> nrow;
        *win_size = new int[2];
        (*win_size)[0] = ncol; //use XSize x YSize instead of nrow x ncol
        (*win_size)[1] = nrow;
      }
      if(*offset_1 == NULL){
        (*offset_1) = new int[2];
        cout << "Please enter the X-offset for image 1: ";
        cin >> (*offset_1)[0];
        cout << "Please enter the Y-offset for image 1: ";
        cin >> (*offset_1)[1];
      }
      if(*offset_2 == NULL){
        (*offset_2) = new int[2];
        cout << "Please enter the X-offset for image 2: ";
        cin >> (*offset_2)[0];
        cout << "Please enter the Y-offset for image 2: ";
        cin >> (*offset_2)[1];
      }
    }

    /*************************************************************************/
    void getOutputFileInfo(string& output_file, string& format){
      if(format.empty()){
        std::cout<< "Valid formats: geotransformiff, PCIDSK, HFA, ENVI" << std::endl;
        std::cout<< "Please enter output format: ";
        getline(std::cin, format);
      }
      if(output_file.empty()){
        std::cout<< "Please enter output file name";
        getline(std::cin, output_file);
      }
    }

    /*************************************************************************/

    //TODO: NEEDS REWRITING!

    // void writeOutputToFile(GDALDataset* outfile, double* tile,
    //                        MatrixXd& A, MatrixXd& B, //Eigenvector matrices
    //                        int xoffset, int yoffset, int ncol, int nrow,
    //                        GDALRasterBand** bands_1,
    //                        GDALRasterBand** bands_2,
    //                        GDALDataset* reference_file,
    //                        VectorXd& sigMADs){
    //   int nBands = A.cols();
    //   double* geotransform   = new double[6];
    //   reference_file->GetGeoTransform(geotransform);
    //     //Move the origins of the picture to the overlap zone. Not implemented yet.
    //       //  geotransform[0] = geotransform[0] + xof*geotransform[1]
    //       //  geotransform[3] = geotransform[3] + y10*geotransform[5]
    //   const char* projection = reference_file->GetProjectionRef();
    //   outfile->SetGeoTransform(geotransform);
    //   outfile->SetProjection(projection);
    //
    //   vector<GDALRasterBand*> outbands  = vector<GDALRasterBand*>(nBands);
    //
    //   for(int i = 0; i < nBands; i++){
    //     outbands[i] = outfile->GetRasterBand(i+1);
    //   }
    //
    //   for(int row = 0; row < nrow; row++){
    //     for(int k = 0; k < nBands; k++){
    //       bands_1[k]->RasterIO(GF_Read, xoffset, yoffset + row, ncol, 1,
    //                            &tile[k], ncol, 1,
    //                            GDT_Float64, sizeof(double)*(nBands*2), 0);
    //       bands_2[k]->RasterIO(GF_Read, xoffset, yoffset + row, ncol, 1,
    //                            &tile[k + nBands], ncol, 1,
    //                            GDT_Float64, sizeof(double)*(nBands*2), 0);
    //     }
    //     MapRMMatrixXd tileMat(tile, ncol, 2*nBands);
    //     MatrixRXd top = tileMat.block(0,0,ncol,nBands);
    //     MatrixRXd bot = tileMat.block(0,nBands,ncol,nBands);
    //     MatrixXd mads = (top * A) - (bot * B);
    //     for(int k = 1; k < nBands; k++){
    //       outbands[k]->RasterIO(GF_Write, xoffset, yoffset + row, ncol, 1,
    //                            &(mads.data()[k*ncol]), ncol, 1,
    //                            GDT_Float64, 0,0 );
    //       //implement chisqr to last band
    //       outbands[k]->FlushCache();
    //     }
    //     imad_utils::rowwise_divide(mads,sigMADs);
    //     //Take columnwise sum of squares, result is (1 x nBands) row vector
    //     VectorXd chisqr = mads.array().square().rowwise().sum().matrix();
    //   }
    //   delete[] geotransform;
    // }
}
