#ifndef HELPERS_H
#define HELPERS_H

#include "gdal_priv.h"
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>


using namespace Eigen;

inline int min(int a, int b){  return a < b ? a : b; }

/* This file contains the headers for the iMad project. */
/* All GDAL functions derive from gdal_priv.h           */
/* Namespaces are named after their .cpp files, e.g. all functions in namespace
 * GdalFileIO are implemented in GdalFileIO.cpp. */

typedef Map<Matrix<double,Dynamic,Dynamic,RowMajor> > MapRMMatrixXd;
typedef Matrix<double,Dynamic,Dynamic,RowMajor> MatrixRXd;

namespace GdalFileIO{

  GDALDataset* openFile(std::string filename);
  int* selectBands(int nBands);
  void fix_missing_band_data(int** bands1, int** bands2, int& nBands);
  void fix_missing_dims_data(int** win_size, int** offset_1,int** offset_2);
  bool has_errors(GDALDataset* file1, GDALDataset* file2,
                 int* bands1_arg, int* bands2_arg, int nBands,
                 int* win_size, int* offsets_1, int* offsets_2, int inp_pen );
  void getOutputFileInfo(std::string& output_file, std::string& format);
  void writeOutputToFile(GDALDataset* outfile, double* tile,
                         MatrixXd& A, MatrixXd& B, //Eigenvector matrices
                         int x10, int y10, int x20, int y20,
                         int ncol, int nrow, int nb_pr, int bufsize,
                         GDALRasterBand** bands_1,
                         GDALRasterBand** bands_2,
                         GDALDataset* reference_file,
                         VectorXd& sigMADs);
}

namespace imad_utils{

  struct Eigentup{
    bool operator < (Eigentup const &other) const;
    Eigentup(double val, VectorXd vec_a, VectorXd vec_b);

    double eigenval;
    VectorXd eigenvec_a;
    VectorXd eigenvec_b;
  };

  void reorder_eigens(VectorXd& lambda, MatrixXd& A, MatrixXd& B);
  void rowwise_subtract(MatrixRXd& A, VectorXd& toSubtract1,
                        MatrixRXd& B, VectorXd& toSubtract2);
  void rowwise_divide(MatrixXd& A, VectorXd& toDivide);
  void colwise_multiply(MatrixXd& A, VectorXd& toMult);

}

namespace geo_utils{

  struct ImageInfo{
    int ncol;
    int nrow;
    int nBands;
    char* projection;
    double* geotransform;
    ImageInfo(GDALDataset* target);
    ~ImageInfo();

    bool operator == (const ImageInfo& other) const;
    bool compatible (const ImageInfo& other) const;
  };

  class CoordTransform{
     Matrix3d Img2Geo;
     Matrix3d Geo2Img;
     Vector3d  input;
     Vector3d  output;

   public:
     //Needed for vectorizable ops
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CoordTransform(double* transform_coeff);
    CoordTransform(GDALDataset* input_file);
    ~CoordTransform();

    //return one coordinate
    double ImgtoGeo_X(double imgX, double imgY);
    double ImgtoGeo_Y(double imgX, double imgY);
    double GeotoImg_X(double geoX, double geoY);
    double GeotoImg_Y(double geoX, double geoY);

    //Modify coordinates in-place
    void GeotoImg(double& imgX, double& imgY);
    void ImgtoGeo(double& geoX, double& geoY);
  };
}

namespace imad_bigfun{

  VectorXd& calc_weights(double* tile, VectorXd& weights, MatrixXd& A,
                        MatrixXd &B, VectorXd& means1, VectorXd& means2,
                        VectorXd& sigMADs, int ncol, int nBands);
  int find_chunksize(int& buffersize, int ncol, int nBands);
  void math_cleanup(MatrixXd& A, MatrixXd& B, MatrixXd& s11, MatrixXd& s12);
void readToBuf(double* tile, GDALRasterBand* band,
        int xstart, int ystart, int bufsize, int nBands);

}

class ImageStats{
  int n2Bands;
  double sum_weights;
  VectorXd means;
  MatrixXd covar;

 public:
  ImageStats(size_t bands);
  ~ImageStats();
  VectorXd get_means();
  MatrixXd get_covar();

  void reset();
  void update(double* input, double* weights, int nrow, int ncol);
};

#endif
