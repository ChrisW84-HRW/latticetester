// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2018  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/*
 * The purpose of this example is to compare the time
 *  of each reduction method. We compute several
 * random matrix and apply algorithms. The parameters of this
 * analysis are defined at the beginning of this file,
 * after include part. The output is a graphic generated by R.
 * The User need R distribution and the header RInside.h
 */


// Include Header
#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

// Include LatticeTester Header
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/Reducer.h"
#include "latticetester/ParamReader.h"
//#include "latticetester/Config.h"
#include "../examples/LatticeAnalysis.h"
//#include "LatticeAnalysis.h"
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkL2.h"
#include "latticetester/NormaMinkL1.h"
#include "latticetester/NormaPalpha.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/Reducer.h"
#include "latticetester/Writer.h"
#include "latticetester/WriterRes.h"
#include "../examples/Config.h"
#include "latticetester/ParamReader.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/NTLWrap.h"


// Include NTL Header
#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "NTL/vec_ZZ.h"
#include "NTL/vec_ZZ_p.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>

using namespace std;
using namespace LatticeTester;

int ZZRR(int argc, char* argv[]) {
  LatticeAnalysis< NTL::ZZ, NTL::RR, NTL::RR> latAnalysis;
  struct stat buf; // properties of a file or directory
  int status = 0;

  for (int j = 2; j < argc; j++) {
    // Do the test for each data file or directory on the command line

    stat(argv[j], &buf);
    if (0 != S_ISDIR(buf.st_mode)) // directory
      status |= latAnalysis.doTestFromDirectory (argv[j]);
    else {
      string dataname(argv[j]);
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);
      status |= latAnalysis.doTestFromInputFile (argv[j]);
    }
  }
  return 1;
}
int ZZDD(int argc, char* argv[]) {
  LatticeAnalysis<NTL::ZZ, double, double> latAnalysis;
  struct stat buf; // properties of a file or directory
  int status = 0;

  for (int j = 2; j < argc; j++) {
    // Do the test for each data file or directory on the command line

    stat(argv[j], &buf);
    if (0 != S_ISDIR(buf.st_mode)) // directory
      status |= latAnalysis.doTestFromDirectory (argv[j]);
    else {
      string dataname(argv[j]);
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);
      status |= latAnalysis.doTestFromInputFile (argv[j]);
    }
  }
  return 1;
}
int LLDD(int argc, char* argv[]) {
  LatticeAnalysis< std::int64_t, double, double> latAnalysis;
  struct stat buf; // properties of a file or directory
  int status = 0;

  for (int j = 2; j < argc; j++) {
    // Do the test for each data file or directory on the command line

    stat(argv[j], &buf);
    if (0 != S_ISDIR(buf.st_mode)) // directory
      status |= latAnalysis.doTestFromDirectory (argv[j]);
    else {
      string dataname(argv[j]);
      dataname.append(".dat");
      stat(dataname.c_str(), &buf);
      status |= latAnalysis.doTestFromInputFile (argv[j]);
    }
  }
  return 1;
}
//==============================================================================

int main (int argc, char *argv[])
{
  std::clock_t start_t = std::clock();
  if (argc < 2) {
    cerr << "\n*** Usage:\n   "
      << argv[0] << " data_file1 data_file2 ...." << endl
      << "or\n   "
      << argv[0] << " dir1 dir2 ...." << endl
      << endl;
    return -1;
  } else if (argc < 3) {
    cerr << "Must specify types\n";
    return -1;
  }
  //cout << "ARGC " << argc << endl;
  if (!strcmp(argv[1], "ZZRR")) {
    ZZRR(argc, argv);
    std::clock_t end_t = std::clock();
    cout << "Time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC << " sec\n";
  } else if (!strcmp(argv[1], "ZZDD")) {
    ZZDD(argc, argv);
    std::clock_t end_t = std::clock();
    cout << "Time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC << " sec\n";
  } else if (!strcmp(argv[1], "LLDD")) {
    LLDD(argc, argv);
    std::clock_t end_t = std::clock();
    cout << "Time: " << (double)(end_t-start_t)/CLOCKS_PER_SEC << " sec\n";
  } else {
    cerr << "No valid type name has been specified\n";
  }

  return 0;
}


