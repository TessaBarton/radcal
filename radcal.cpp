/*H**********************************************************************
* FILENAME :        radcal.cpp
*
* DESCRIPTION :
*       Takes in matrix of weights as input. Determines top n% of weights to
b*e used as points of no change for radiometric normalization. Using the (x,y
) coordinates of these
*
* PUBLIC FUNCTIONS :
*       type    FunctionName( Parameter1, Parameter2 )
*
* NOTES :
*
*
*       Copyright (c) Mort Cantey(2007), Theresa Barton(2014).
*
* CHANGES :
*
*
*H*/

/**************************************************************************
*   INCLUDE FILES
***************************************************************************/

// to compile g++ -std=c++11 -I /usr/local/Cellar/eigen/3.2.1/include/eigen3 radcal.cpp GdalFileIO.cpp -lgdal


#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "stdio.h"
#include "time.h"
#include "radcal.h"
#include "tgmath.h"
// header files from Donovan Parks
#include "Precompiled.hpp"
#include "LinearRegression.hpp"

#include <ogr_spatialref.h>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include <algorithm>



#include <cmath>
#include <algorithm>

/**************************************************************************
*   NameSpaceIng
***************************************************************************/

using namespace std;
using Eigen::MatrixXd;//this represents a matrix of arbitrary size (hence the X in MatrixXd), in which every entry is a doubl
using Eigen::VectorXd;




/**************************************************************************
*   Data points input
***************************************************************************/


/**************************************************************************
*   Linear Regression Code (from Donovan Parks)
***************************************************************************/

bool LinearRegression::LeastSquaresEstimate(std::vector<double>& x,
std::vector<double>& y, RESULTS& results)// needs to return vector of gain and offset information
{
   int numPts = x.size();
   if(numPts < 2)
      return false;

      double sumX = 0;
      double sumY = 0;
      double sumXY = 0;
      double sumX2 = 0;
      double sumY2 = 0;
      for(int i=0; i < numPts; ++i)
         {
            sumX += x.at(i);
            sumY += y.at(i);
            sumXY	+= x.at(i) * y.at(i);
            sumX2	+= x.at(i) * x.at(i);
            sumY2	+= y.at(i) * y.at(i);
         }

         double SSxy = sumXY - (sumX*sumY) / numPts;
         double SSxx = sumX2 - (sumX*sumX) / numPts;
         double SSyy = sumY2 - (sumY*sumY) / numPts;

         if(SSxx != 0 && SSyy != 0)
            {
               results.slope = SSxy / SSxx;

               results.cod = (SSxy*SSxy) / (SSxx*SSyy);
               results.r = SSxy / sqrt(SSxx*SSyy);
            }
            else if(SSxx == 0 && SSyy == 0)
            {
               results.slope = DBL_MAX;
               results.cod = DBL_MAX;
               results.r = 1;
            }
         else
            {
               results.slope = DBL_MAX;

               results.cod = 1;
               results.r = -1;
            }

            results.offset = (sumY - results.slope*sumX) / numPts;
            results.dof = numPts-2;
            results.nPts = numPts;

            double SSE = SSyy-results.slope*SSxy;
            results.sd = sqrt(SSE / results.dof);
            results.seos = results.sd / sqrt(SSxx);
            results.seoi = results.sd*sqrt((1/numPts) + ((sumX/numPts)*(sumX/numPts))/SSxx);

            results.tStat = results.slope / results.seos;

            return true;
         }

         /**
         * Calculate the error of a given best fit line.
         */
         double LinearRegression::Error(int numPts, double sumX, double sumY, double sumXY,
         double sumX2, double sumY2)
         {
            double SSxx = sumX2 - (sumX*sumX) / numPts;
            double SSyy = sumY2 - (sumY*sumY) / numPts;
            double SSxy = sumXY - (sumX*sumY) / numPts;

            double slope = SSxy / SSxx;
            double offset = (sumY - slope*sumX) / numPts;

            return SSyy - slope*SSxy;
         }

         /**************************************************************************
         *   image alterations
         ***************************************************************************/
         void normalize(GDALDataset* file,MatrixXd gainsandoffsets){

           //vars
           int MaxVal = 255;
           string format      = "GTiff";
           string output_file = "output_file";
           int ncol        = file->GetRasterYSize();
           int nrow        = file->GetRasterXSize();
           int nBands      = file->GetRasterCount();

           //create outfile
           GdalFileIO::getOutputFileInfo(output_file, format);
           GDALDriver* outdriver =GetGDALDriverManager()->GetDriverByName(format.c_str());
           GDALDataset *outfile  = outdriver->Create(output_file.c_str(),ncol,nrow,
           nBands, GDT_Float64, NULL);
           int numbands= file->GetRasterCount();
           GDALRasterBand *inputband;
           GDALRasterBand *outputband;

           //create buffer
           double * tile= new double[ncol];

           // Main loop
           for( int band=1;band<numbands+1;band++){
             inputband = file->GetRasterBand(band);
             outputband = outfile->GetRasterBand(band);
             cout << "fetched bands\n";
             double gain   = gainsandoffsets(band-1,0);
             double offset = gainsandoffsets(band-1,1);
             for(int row = 0;row<nrow; row++){
               inputband->RasterIO( GF_Read, 0, row, ncol, 1,tile, ncol, 1, GDT_Float64,0,
               0);
               //cout << "got input band =>";

               for(int pixel =0; pixel < ncol; pixel ++){
                 double value = tile[pixel];
                 //cout << value;
                 if (value>1){ // so 1 means no data need to handle this!
                   value = value*gain + offset;
                   if(value>MaxVal){
                     value = MaxVal;
                   }
                   tile[pixel]=value;
                 }
                   else{
                     tile[pixel]=1;
                   }
                   //cout << value << "\n";
                 }

                outputband->RasterIO( GF_Write, 0, row, ncol, 1,tile, ncol, 1, GDT_Float64,0,
                   0);
                   //cout << "wrote output band \n";
             }
           }

           delete[] tile;
         }



/**************************************************************************
*   Data points input
***************************************************************************/
struct inputPixel {
  double radiance;
  int x;
  int y;
};
vector<pair<int,int> > * xyget(string& filename,double threshold){ // nielson paper 2003
vector<pair<int,int> > * generatedpoints = new vector<pair<int,int> >();
   generatedpoints->push_back(make_pair(1,1));
   generatedpoints->push_back(make_pair(9733	,11170	));
   generatedpoints->push_back(make_pair(10946, 3574));

   generatedpoints->push_back(make_pair(13862,12724));
   generatedpoints->push_back(make_pair(6520,10070));
   generatedpoints->push_back(make_pair(16168, 13875));
   generatedpoints->push_back(make_pair(8449, 5886));
   generatedpoints->push_back(make_pair(479, 6309));
   generatedpoints->push_back(make_pair(1213,5256));
   generatedpoints->push_back(make_pair(10283,1856));
   return generatedpoints;
//   vector<pair<int,int> > * staticPoints = new vector<pair<int,int> >();
//
//   // open file
//   string maskImage   = "C++_asd.tif";
//   GDALDataset* file  = GdalFileIO::openFile(maskImage);
//   int ncol        = file->GetRasterXSize();
//   int nrow        = file->GetRasterYSize();
//   int lastband    = file->GetRasterCount();
//   //cout << lastband;
//   int histogram[5];
//
//   /*initialize array of structs, array of doubles*/
//   inputPixel * inputImage = new inputPixel[nrow*ncol];
//
//   GDALRasterBand *inputband= file->GetRasterBand(lastband);
//   double * tile= new double[ncol];
//   // read in change probabilities
//   for(int row = 0;row<nrow; row++){
//     inputband->RasterIO( GF_Read, 0, row, ncol, 1,tile, ncol,
//      1, GDT_Float64,0,0);
//
//     for(int col= 0; col<ncol;col++){
//         // cout << tile[col]<< ",";
//         // if (tile[col] >=0 && tile[col] <10) {histogram[0]++; }
//         // else if (tile[col] >=10 && tile[col] <100) {histogram[1]++; }
//         // else if (tile[col] >= 100 && tile[col] <1000) {histogram[2]++; }
//         // else if (tile[col] >=100 && tile[col] <1000) {histogram[3]++; }
//         // else if (tile[col] >=1000 && tile[col] <10000) {histogram[4]++; }
//         // else if (tile[col] >=10000 && tile[col] <100000) {histogram[5]++; }
//       inputImage[row + col].radiance = tile[col];
//       inputImage[row+col].x = col;
//       inputImage[row+col].y = row;
//     }
//     }
//     // sort probabilites
//     sort(inputImage,inputImage +(nrow*ncol), [](const inputPixel &a, const inputPixel  &b){ return a.radiance < b.radiance; });
//
//     // //histogram of values, hmm chi squared interface?
//     // cout<< "0:10 ="<< histogram[0];
//     // //for (int i=0;i<histogram[0];i++) {cout << "*"; }
//     // cout << "\n"<< "10:100 = "<< histogram[1];
//     // //for (int i=0;i<histogram[1];i++) {cout << "*"; }
//     // cout << "\n"<< "100:1000 = "<< histogram[2];
//     // //for (int i=0;i<histogram[2];i++) {cout << "*"; }
//     // cout << "\n"<< "1000:10000 = "<< histogram[3];
//     // //for (int i=0;i<histogram[3];i++) {cout << "*"; }
//     // cout << "\n"<< "10000:100000 = "<< histogram[4]<< "\n";
//     // //for (int i=0;i<histogram[4];i++) {cout << "*"; }
//
//     // cout << inputImage[0].radiance << ",";
//     // cout << inputImage[ncol*nrow].radiance;
//
// // still need to cut off values that are below threshold
// // still need to put the xy coords of those values in static points
//
//     delete[] tile;
//     delete[] inputImage;
//
//    return staticPoints;
 }



         /**************************************************************************
         *   Radcal
         ***************************************************************************/
         void radcal(){


            GDALAllRegister();
            OGRSpatialReference oSRS;
            // inputs
            string filename="tjpeg.tif";
            string filenameagain = "tjpeg.tif";
            //extract bands at (x,y)
            GDALDataset* file1 = GdalFileIO::openFile(filename);
            GDALDataset* file2 = GdalFileIO::openFile(filenameagain);


            // find xy of the points above threshold t
            double threshold = .05;
            string maskImage = "C++_asd.tif";
            vector<pair<int,int> > *generatedpoints = xyget(maskImage,threshold);// refactor to static points





            int   nXSize = file1->GetRasterXSize();
            int   nYSize = file1-> GetRasterYSize();
            int numpoints =generatedpoints->size();
            int numbands = file1->GetRasterCount();
            MatrixXd pixValsMat(numbands,numpoints); // need to change to doubles
            GDALRasterBand *rasterband;
            int nBlockXSize, nBlockYSize;

         // read in data from image one at points

            for( int band=1;band<numbands+1;band++){
               //cout << "_____" << "band" << i << "____"<<"\n";
               for(int point = 0;point<numpoints; point++){
                  rasterband = file1->GetRasterBand(band);
                  int pixel= generatedpoints-> at(point).first;
                  int line = generatedpoints-> at(point).second;
                  int nXSize = rasterband->GetXSize();
                  double * value;

                  value = (double *) CPLMalloc(sizeof(double));

                  rasterband->RasterIO( GF_Read, pixel, line, 1, 1,value, 1 /*nXSize*/, 1, GDT_Float64,
                  0, 0 );
                  pixValsMat(band-1, point)=*value;
                  //cout << *value<< "\n";

                  CPLFree(value);
               }
            }

            // read in data from image 2 at points

            MatrixXd pixValsMat2(numbands,numpoints);
            for( int band=1;band<numbands+1;band++){
               //cout << "_____" << "band" << i << "____"<<"\n";
               for(int point = 0;point<numpoints; point++){
                  rasterband = file2->GetRasterBand(band);
                  int pixel= generatedpoints-> at(point).first;
                  int line = generatedpoints-> at(point).second;
                  int nXSize = rasterband->GetXSize();
                  double * value;
<<<<<<< HEAD
                  value = (double *) CPLMalloc(sizeof(double));
=======
                  value = (double *) CPLMalloc(sizeof(double)*nXSize);
>>>>>>> c76a2baae3282ae5aa52f7de10d84622c052f8a7
                  rasterband->RasterIO( GF_Read, pixel, line, 1, 1,value, 1 /*nXSize*/, 1, GDT_Float64,
                  0, 0 );
                  pixValsMat2(band-1, point)=*value;
                  //cout << *value<< "\n";

                  CPLFree(value);
               }
            }


            /*run Regression*/
            MatrixXd lrparams(numbands,2);

            // vectors are rows, pass vectors to linreg
            for(int band =0; band < numbands; band++){

               LinearRegression lin =LinearRegression();
               double pixValsArr [numpoints];
               for( int pixel = 0; pixel <numpoints; pixel++){
                  pixValsArr[pixel] = pixValsMat(band,pixel);
               }
               vector<double> pixValsVec (pixValsArr, pixValsArr + sizeof(pixValsArr) / sizeof(double) );

               double pixValsArr2[numpoints];
               for( int pixel = 0; pixel <numpoints; pixel++){
                  pixValsArr2[pixel] = pixValsMat2(band,pixel);
               }
               vector<double> pixValsVec2 (pixValsArr2, pixValsArr2 + sizeof(pixValsArr2) / sizeof(double) );
               LinearRegression::RESULTS results = LinearRegression::RESULTS();
               lin.LeastSquaresEstimate(pixValsVec, pixValsVec2,results);//member function
               double gain = results.slope;
               cout << "gain for band n : " <<gain <<"\n";
               double offset = results.offset;
               cout << "offset for band n : " <<offset <<"\n";
               lrparams(band,0)=gain;
               lrparams(band,1)=offset;

            }
         for(int x=0;x<numbands;x++){
               for(int y=0;y<2;y++)  // loop for the three elements on the line
                  {
                     cout<<lrparams(x,y)<<" ";  // display the current element out of the array
                  }
                  cout<<endl;  // when the inner loop is done, go to a new line
               }

         // normalize file 2 using params
         normalize(file2,lrparams);
         GDALClose(file1);
         GDALClose(file2);

         }

               int main()// probably should take stuff like threshhold value, image 1, image 2 etc.
               {

                  radcal();
                  return 0;
               }
