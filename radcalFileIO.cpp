//GDALIOFile
#include "gdal_priv.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>

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

    bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2){
      /* Check that dimensions of the two files match in X, Y, and num. bands */
      if(dataset1->GetRasterXSize() != dataset2->GetRasterXSize()){
        cout << "X-dimension mismatch in images. Exiting..." << endl;
        return false;
      }
      if(dataset1->GetRasterYSize() != dataset2->GetRasterYSize()){
        cout << "Y-dimension mismatch in images. Exiting..." << endl;
        return false;
      }
      if(dataset1->GetRasterCount() != dataset2->GetRasterCount()){
        cout << "Raster count mismatch in images. Exiting..." << endl;
        return false;
      }
      //All dimensions match!
      return true;
    }

    vector<int>* selectBands(){
      vector<int>* bands = new vector<int>();
      cout << "Please enter the band numbers you would like to use, followed by a return"
           << endl << "When you are finished, please enter 0. Yes this interface sucks." << endl;

      string input;
      cin >> input;
      while(input != "0"){
        bands->push_back(atof(input.c_str()));
        cin >> input;
      }
      if (bands->size() < 1){
        throw std::out_of_range("Need to specify at least one band.");
      }
      return bands;
    }

}