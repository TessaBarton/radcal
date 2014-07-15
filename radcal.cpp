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

#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXf;//this represents a matrix of arbitrary size (hence the X in MatrixXd), in which every entry is a doubl
using Eigen::VectorXf;

vector<pair<int,int> > * xygen(int n){
   vector<pair<int,int> > * points = new vector<pair<int,int> >();
   for(int i =0; i< n; i++){
    points->push_back(make_pair(i,i));
    //cout << points-> at(i).first << ", " << points-> at(i).second << '\n';
   }
   return points;

}


void radcal(){
    
	GDALAllRegister();
	// inputs
	    string filename="tjpeg.tif";
	    vector<pair<int,int> > * generatedpoints = xygen(10);
		int numpoints = generatedpoints->size();

	//extract bands at (x,y)
    	GDALDataset* file1 = GdalFileIO::openFile(filename);
    	int numbands = file1->GetRasterCount();
    	cout << numbands;
    
    //for each point read in bands
    

}

int main()
{
  radcal();
  return 0;
}
