//header file
#ifndef HELPERS_H
#define HELPERS_H

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>

/* This file contains the headers for the iMad project. */
/* All GDAL functions derive from gdal_priv.h           */

namespace GdalFileIO{
      void getOutputFileInfo(std::string& output_file, std::string& format);
  GDALDataset* openFile(std::string filename);
  bool dimensionsmatch(GDALDataset* dataset1, GDALDataset* dataset2);
  std::vector<int>* selectBands();

}



#endif
