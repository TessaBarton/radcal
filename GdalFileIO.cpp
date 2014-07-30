
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




}
