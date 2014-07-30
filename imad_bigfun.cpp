#include "gdal_priv.h"
#include "imad.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions.hpp> //for cdf()

using namespace Eigen;

namespace imad_bigfun{

  /* Uses the info given to read an appropriate section of the image into
   * the buffer tile. tile must be at least bufsize large. The chunk returned
   * is a single row at (xoffset + bufsize * nbuf_so_far, yoffset + row) with
   * xsize of bufsize */

  void readToBuf(double* tile, GDALRasterBand* band,
          int xstart, int ystart, int bufsize, int nBands){

    band->RasterIO(GF_Read, xstart, ystart, bufsize, 1,
                       tile, bufsize, 1,
                       GDT_Float64, sizeof(double)*(nBands*2), 0);

  }

  /*************************************************************************/


  /* Gets CDF of Chisquared and places into weights array. CDF function is
   * backwards compared to Python (in python, param to chi2.cdf is second arg)
   * so [C++] cdf(rng(a), b) = stats.chi2.cdf(b,a) [Python] */

   VectorXd& calc_weights(double* tile, VectorXd& weights, MatrixXd& A,
                         MatrixXd &B, VectorXd& means1, VectorXd& means2,
                         VectorXd& sigMADs, int bufsize, int nBands){
     MapRMMatrixXd tileMat(tile, bufsize, 2*nBands);
     //Separate data from images 1 and 2 into their own matrices
     MatrixRXd top = tileMat.block(0,0,bufsize,nBands);
     MatrixRXd bot = tileMat.block(0,nBands,bufsize,nBands);

     imad_utils::rowwise_subtract(top, means1, bot, means2);
     //Multiply images by eigenvectors to get MAD variates
     MatrixXd mads = (top * A) - (bot * B);

     imad_utils::rowwise_divide(mads,sigMADs);

    //Do a rowwise sum of squares to get the chisqr values
     VectorXd chisqr = mads.array().square().rowwise().sum().matrix();

     boost::math::chi_squared dist(nBands);
     for(int i = 0; i < weights.rows(); i++){
       weights(i) = 1 - boost::math::cdf(dist, chisqr(i));
     }

     return weights;
   }

  /*************************************************************************/

   /* Modifies the buffersize in-place, return value is number of chunks n
   needed in order to competely cover the input row. Can theoretically
   handle an image the size of the solar system...hopefully that won't be
   needed anytime soon. */

  int find_chunksize(int& buffersize, int ncol, int nBands){
    const int max_memory_in_bytes = 1000000000; //1GB
    const int max_bufsize = max_memory_in_bytes / ncol / nBands / sizeof(double);

    if(ncol < max_bufsize){
      buffersize = ncol;
      return 1;
    }
    else{
      for(int i = 1; i < std::numeric_limits<int>::max() - 2; i++){
        if(ncol / i < max_bufsize){
          buffersize = ncol / i + 1; //rounding up
          return i;
        }
      }
    }
    throw std::invalid_argument("find_chunksize: Image is too large to chunk...what the hell do you have?");
  }

  /*************************************************************************/

  /* What follows are some of the most poorly-commented lines from the original
  * Python implementation, and my understanding of them is still fuzzy.
  * I have included modified versions of the original comments, and the
  * lines are each numbered, with comment-footnotes at the end of the file */

  //You are not expected to understand these lines. I certainly don't.

  void math_cleanup(MatrixXd& A, MatrixXd& B, MatrixXd& s11, MatrixXd& s12){
  /* 1 */ MatrixXd D = s11.diagonal().array()      //Bookkeeping
                        .sqrt().cwiseInverse()   //Inverse sqrt
                        .matrix().asDiagonal();  //More bookkeeping
  /* 2 */ VectorXd s = (D*s11*A).rowwise().sum();
        VectorXd tmp1 = (s.array() / (s.array().abs())).matrix();
  /* 3 */ imad_utils::colwise_multiply(A, tmp1);
  /* 4 */ VectorXd cov = (A.transpose() * s12 * B).diagonal();
  /* 5 */ MatrixXd tmp2 = (cov.array() / cov.array().abs()).matrix().asDiagonal();
        B = B * tmp2;
  }

//End namespace
}

/**** Notes on tile matrix ****/
/* tile is a strange name for a matrix (borrowed from the python code).
 * It holds a single row of both images. The top half of the matrix holds
 * the row from file1, with all of its bands in separate rows (e.g. for any
 * given row, row 1 of tile will hold image1, band1, row 2 will hold image1,
 * band2, etc.). The same is true of the second image in the bottom
 * half of tile. Note that this tile is transposed compared to the one in
 * Python: Pyversion has rows in cols and files split left/right. */

/**** Notes on numbered lines in imad() ****/
/* While going through these, it is helpful to remember that left multiplying by
 * a diagonal matrix (D*A) multiplies each row by the diagonal element, i.e.
 * row1(A) = row1(A) * D(1,1), row2(A) = row2(A) * D(2,2), etc. and right-multiplying
 * by a diagonal matrix multiplies the columns instead of the rows. Off we go!*/

/* 1 */ /* Completely unintelligible. First, we take the diagonal of s11. We then
 * cast it to an array (to deal with some issues with eigen), then take an
 * element-wise inverse square root. We then revert it back to a diagonal matrix.
 * The end effect is identical to masking all but the diagonal of s11, then
 * taking an element-wise inverse square root. */

/* 2 */ /* s is actually a row-vector. Not sure what this actually does */

/* 3 */ /* If S(i,i) was negative, we set it to -1. Otherwise, we set it to +1.
 * We then multiply A by this matrix to ensure the correlation sum is positive
 * (Yes, this part feels like magic to me too). */

/* 4,5 */ /* A similar process for B */
