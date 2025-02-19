// File `testMatrixCreationSpeed`

/**
 * This program compares the effectiveness of different resizing options when using
 * an `IntMat` object whose required size changes frequently.
 * One option is to create a single matrix of sufficiently large dimension for all
 * cases and never resize it, but just use the upper left part for the size that we need.
 * Other options are to resize only the number of rows, or resize the full matrix only when
 * the required dimension changes if this is not too frequent,
 * or to resize the full matrix each time, or to create a new matrix object
 * each time. We compare all these methods, in dimensions that vary from 5 to 50.
 */

#define TYPES_CODE  LD     // Int = int64_t
// #define TYPES_CODE  ZD     // Int = ZZ

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"

using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are will be passed as template parameters from the `main`.
const int64_t maxNumSizes = 6; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 20, 30, 40, 50 };
const int64_t numMeth = 5; // Number of different methods to be tested.
clock_t timer[numMeth][maxNumSizes+1]; // timer[meth][0] used for total over dims in swapLoops version.
Int sum[numMeth][maxNumSizes+1];
std::string names[numMeth] = {
      "One matrix, no resize         ",
      "One matrix, resize rows often ",
      "One matrix, resize when needed",
      "One matrix, resize all often  ",
      "New matrix for each repet.    " };

void printTables(int64_t numSizes) {
   int64_t d, total;
   std::cout << "Dimension:                     ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "  ";
   std::cout << "    total \n\n";
   std::cout << "Timings for the different tasks, in basic clock units (microseconds): \n";
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      total = 0;
      for (d = 0; d < numSizes; d++) {
         std::cout << std::setw(9) << timer[meth][d] << " ";
         total += timer[meth][d];
      }
      std::cout << std::setw(9) << total << " ";
      std::cout << "\n";
   }
   std::cout << "\n";
   std::cout << "Sums of diagonal entries:\n";
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      for (d = 0; d < numSizes; d++)
         std::cout << std::setw(9) << sum[meth][d] << " ";
      std::cout << "\n";
   }
   std::cout << "\n";
}

void printTablesSwapLoops() {
   std::cout << "Timings (total) in basic clock units (microseconds): \n";
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      std::cout << std::setw(9) << timer[meth][0] << "\n";
   }
   std::cout << "\n";
}

// Testing loop. For each matrix size, we fill the diagonal `numReps' times.  We do that in 5 different ways.
void testLoop(const int64_t &numReps, const int64_t &numSizes) {
   clock_t tmp;
   int64_t i, j, d, dim;
   int64_t maxdim = dimensions[numSizes - 1];  // The largest dimension that we actually use.
   IntMat A, B;

   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes   
   std::cout << "**************************************************************\n";
   std::cout << "TestMatrixCreationSpeed, outside loop is on dimension, inside loop on replications. \n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of repetitions: " << numReps << "\n";

   // Method 1: Create a single matrix object of maximal dimension, never resize.
   A.SetDims(maxdim, maxdim);
   for (d = 0; d < numSizes; d++) {
      dim = dimensions[d];
      tmp = clock();
      for (i = 0; i < numReps; i++) {
         for (j = 0; j < dim; j++)
            A[j][j] = i + j;
      }
      timer[0][d] = clock() - tmp;
      sum[0][d] = 0;
      for (j = 0; j < dim; j++)  // This summation part is not timed.
         sum[0][d] += A[j][j];
   }

   // Method 2: Create a single matrix object, and resize only the number of rows, each time.
   for (d = 0; d < numSizes; d++) {
      dim = dimensions[d];
      tmp = clock();
      for (i = 0; i < numReps; i++) {
         B.SetDims(0, maxdim);
         B.SetDims(dim, maxdim);
         for (j = 0; j < dim; j++)
            B[j][j] = i + j;
      }
      timer[1][d] = clock() - tmp;
      sum[1][d] = 0;
      for (j = 0; j < dim; j++)
         sum[1][d] += B[j][j];
   }

   // Method 3: Create a single matrix object, resize each time, but the resize
   // is effective only when the dimension changes.
   for (d = 0; d < numSizes; d++) {
      dim = dimensions[d];
      tmp = clock();
      for (i = 0; i < numReps; i++) {
         B.SetDims(dim, dim);   // This does something only when the dim has changed.
         for (j = 0; j < dim; j++)
            B[j][j] = i + j;
      }
      timer[2][d] = clock() - tmp;
      sum[2][d] = 0;
      for (j = 0; j < dim; j++)
         sum[2][d] += B[j][j];
   }

   // Method 4: Create a single matrix object, but resize it completely each time.
   for (d = 0; d < numSizes; d++) {
      dim = dimensions[d];
      tmp = clock();
      for (i = 0; i < numReps; i++) {
         B.SetDims(0, 0);
         B.SetDims(dim, dim);
         for (j = 0; j < dim; j++)
            B[j][j] = i + j;
      }
      timer[3][d] = clock() - tmp;
      sum[3][d] = 0;
      for (j = 0; j < dim; j++)
         sum[3][d] += B[j][j];
   }

   // Method 5: Declare and create a new matrix object each time.
   for (d = 0; d < numSizes; d++) {
      dim = dimensions[d];
      tmp = clock();
      for (i = 0; i < numReps; i++) {
         IntMat C;
         C.SetDims(dim, dim);
         sum[4][d] = 0;
         for (j = 0; j < dim; j++) {
            C[j][j] = i + j;
            sum[4][d] += C[j][j];  // Removing this statement does not impact the speed much.
         }
      }
      timer[4][d] = clock() - tmp;
   }

   printTables(numSizes);
}


// Testing loop. For each matrix size, we fill the diagonal `numReps' times.
// Here we exchange the order of the loops.
void testSwapLoops(const int64_t &numReps, const int64_t &numSizes) {
   clock_t tmp;
   int64_t i, j, d, dim;
   int64_t maxdim = dimensions[numSizes - 1];
   IntMat A, B;

   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "**************************************************************\n";
   std::cout << "TestMatrixCreationSpeed, inside loop is on dimension. \n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of repetitions: " << numReps << "\n";

   // Method 1: Create a single matrix object of maximal dimension, never resize.
   tmp = clock();
   A.SetDims(maxdim, maxdim);
   for (i = 0; i < numReps; i++) {
      for (d = 0; d < numSizes; d++) {
         dim = dimensions[d];
         for (j = 0; j < dim; j++) {
            A[j][j] = i + j;
         }
      }
   }
   timer[0][0] = clock() - tmp;

   // Method 2: Create a single matrix object, and resize only the number of rows, each time.
   tmp = clock();
   for (i = 0; i < numReps; i++) {
      for (d = 0; d < numSizes; d++) {
         dim = dimensions[d];
         // B.SetDims(0, maxdim);
         B.SetDims(dim, maxdim);
         for (j = 0; j < dim; j++) {
            B[j][j] = i + j;
         }
      }
   }
   timer[1][0] = clock() - tmp;

   // Method 3: Create a single matrix object, resize each time, but the resize
   // is effective only when the dimension changes.
   tmp = clock();
   for (i = 0; i < numReps; i++) {
      for (d = 0; d < numSizes; d++) {
         dim = dimensions[d];
         B.SetDims(dim, dim);   // This does something only when the dim has changed.
         for (j = 0; j < dim; j++) {
            B[j][j] = i + j;
         }
      }
   }
   timer[2][0] = clock() - tmp;

   // Method 4: Create a single matrix object, but resize it completely each time.
   tmp = clock();
   for (i = 0; i < numReps; i++) {
      for (d = 0; d < numSizes; d++) {
         dim = dimensions[d];
         B.SetDims(0, 0);
         B.SetDims(dim, dim);
         for (j = 0; j < dim; j++) {
            B[j][j] = i + j;
         }
      }
   }
   timer[3][0] = clock() - tmp;

   // Method 5: Declare and create a new matrix object each time.
   tmp = clock();
   for (i = 0; i < numReps; i++) {
      for (d = 0; d < numSizes; d++) {
         dim = dimensions[d];
         IntMat C;
         C.SetDims(dim, dim);
         for (j = 0; j < dim; j++) {
            C[j][j] = i + j;
         }
      }
   }
   timer[4][0] = clock() - tmp;

   printTablesSwapLoops();
}


int main() {

   const int64_t numSizes = 6;
   const int64_t numReps = 100000;

   testLoop(numReps, numSizes);
   testSwapLoops(numReps, numSizes);
   return 0;
}

