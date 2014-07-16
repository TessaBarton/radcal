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
// to compile g++ -I /usr/local/Cellar/eigen/3.2.1/include/eigen3 radcal.cpp GdalFileIO.cpp -lgdal
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "stdio.h"
#include "time.h"
#include "radcal.h"
#include "tgmath.h"

#include <ogr_spatialref.h>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>


using namespace std;
using Eigen::MatrixXf;//this represents a matrix of arbitrary size (hence the X in MatrixXd), in which every entry is a doubl
using Eigen::VectorXf;

vector<pair<int,int> > * xygen(int n, double * padfTransform, int rasterx, int rastery){

   vector<pair<int,int> > * points = new vector<pair<int,int> >();
   double xwest, xeast, ynorth, ysouth;
   xwest  = padfTransform[0] + 0*padfTransform[1] + 0*padfTransform[2];
   xeast = padfTransform[0] + rasterx *padfTransform[1] + 0*padfTransform[2];
   ynorth = padfTransform[3] + 0*padfTransform[4] + 0*padfTransform[5];
   ysouth = padfTransform[3] + 0*padfTransform[4] + rastery*padfTransform[5];

   for(int i =0; i< n; i++){
      int  first, second;
      first = static_cast<int>(xwest + pow(10,i));
      second = static_cast<int>(ynorth +pow(10,i));
      points->push_back(make_pair(first,second));
      //cout << points-> at(i).first << ", " << points-> at(i).second << '\n';
   }
   return points;

}


void radcal(){

   GDALAllRegister();
   OGRSpatialReference oSRS;
   // inputs
   string filename="tjpeg.tif";


   //extract bands at (x,y)
   GDALDataset* file1 = GdalFileIO::openFile(filename);

   double * transform = new double[6];
   file1->GetGeoTransform(transform);
   //cout << transform[0] << '\n' << transform[3]<< '\n'<< transform[1] << '\n' << transform[5]<< '\n';

   vector<pair<int,int> > * generatedpoints = xygen(10,transform,file1->GetRasterXSize(),file1->GetRasterYSize());
   int numpoints = generatedpoints->size();

//extract points
   int numbands = file1->GetRasterCount();
   // put in own func?
   MatrixXf pixVals(numpoints,numbands);

   for( int i =1;i< numbands;i++){
      GDALRasterBand *band;
      band = file1->GetRasterBand(i);
      for(int j = 0;j<numpoints; j++){
         int pixel= generatedpoints-> at(j).first;
         int line = generatedpoints-> at(j).second;
         float val;
         float *buffer;
         file1->RasterIO(GF_Read,pixel,line,0,0,buffer, 1, 1,GDT_Float32,1,NULL, 0, 0,0);
         cout << *buffer;



      }
   }

   // Band->RasterIO(hband, GF_Read, ;


   //for each point read in bands
   delete[] transform;

}

int main()
{
   radcal();
   return 0;
}
