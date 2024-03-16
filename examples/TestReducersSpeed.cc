/**
 * In this example, we compare different ways of computing the shortest vector
 * in a lattice or in its $m$-dual.  For this, we do the BB after some type of
 * pre-reduction of the basis. We compare the total times to do that
 * (pre-red. + BB) for various types of pre-reductions, with lattices that come
 * from Korobov lattice rules in 10 to 40 dimensions, with prime modulus `m`.
 * For the BB, we use the Cholesky decomposition.
 * The pre-reductions considered are:
 * none, pairwise, LLL with delta = 0.5, 0.8, and 0.99999, and BKZ with the
 * default parameters. We do this for two different prime values of `m`.
 * The timings are in terms of total number of clock ticks used by each method.
 *
 * Example of results:
*/

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Select the flexible types Int and Real here.
//#define TYPES_CODE  LD     // Int == int64_t
#define TYPES_CODE  ZD     // Int == ZZ, Real = double
//#define TYPES_CODE  ZR     // Int == ZZ, Real = RR

#include <iostream>
#include <cstdint>
#include <ctime>
#include <type_traits>
#include <typeinfo>

#include <NTL/tools.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"    // This defines Int and Real
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Chrono.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/Rank1Lattice.h"
// #include "latticetester/BasisConstruction.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"

using namespace LatticeTester;

//Int m(1021);     // Modulus m = 1021
Int m(1048573);  // Prime modulus near 2^{20}
//Int m(1073741827);  // Prime modulus near 2^{30}
//Int m(1099511627791);  // Prime modulus near 2^{40}
//Int m(1125899906842597);  // Prime modulus near 2^{50}

const long numSizes = 4;    // Number of matrix sizes (choices of dimension).
const long dimensions[numSizes] = { 10, 20, 30, 40 };
long maxdim = dimensions[numSizes - 1];   // Maximum dimension
const long numMeth = 9;    // Number of methods, and their names.
std::string names[numMeth] = { "LLL5        ", "LLL99999    ",  "BKZ99999    ",
        "Direct BB   ", "Pairwise+BB ", "LLL5+BB     ", "LLL8+BB     ", 
        "LLL99999+BB ", "BKZ99999+BB "};
// PrecisionType prec = QUADRUPLE;  // DOUBLE, QUADRUPLE, XDOUBLE, RR
PrecisionType prec = DOUBLE;

// We use ctime directly for the timings, to minimize overhead.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes]; // Clock times in microseconds.
Real sumSq[numMeth][numSizes];    // Sum of squares of vector lengths (for checking).
RealVec sqlen; // Memory for this vector is reserved in the main.
Rank1Lattice<Int, Real> *korlat;    // We create a single lattice object.
ReducerBB<Int, Real> *red;          // Also a single ReducerBB object.


inline void beforeReduct(long dim) {
    korlat->buildDualBasis(dim);
    korlat->dualize();
    sqlen[0] = 0.0;
    // CopyPartMat (basis2, korlat->getBasis(), dim, dim);  // Copy basis1 to basis2.
    tmp = clock();
    }
    
inline void afterReduct(long meth, long d) {
    timer[meth][d] += clock() - tmp;
    sumSq[meth][d] += sqlen[0];
    std::cout << "Done with method = " << meth << ",  dim = " << dimensions[d] << "\n";
    }   

// Speed test for dim = dimensions[d], with given matrices.
// We test several methods to reduce basis1 and find a shortest vector in its lattice.
static void reduceBasis (long d) {
    long dim = dimensions[d];

    beforeReduct (dim);  // Only pre-reductions.
    redLLLNTL (korlat->getBasis(), &sqlen, 0.5, dim, prec);
    afterReduct (0, d);
 
    beforeReduct (dim);
    redLLLNTL (korlat->getBasis(), &sqlen, 0.99999, dim, prec);
    afterReduct (1, d);

    beforeReduct (dim);
    redBKZ (korlat->getBasis(), &sqlen, 0.99999, 10, 0, dim, prec);
    afterReduct (2, d);
return;

    if (dim < 30) {
      beforeReduct (dim);
      if (!red->shortestVector())
         std::cout << " shortestVector failed with no pre-reduction, dim  = " << dim << "\n";
      afterReduct (3, d);
    }

    if (dim < 30) {
      beforeReduct (dim);
      red->redDieter(0);
      if (!red->shortestVector())
          std::cout << " shortestVector failed for pairwise pre-reduction with dim  = " << dim << "\n";
      afterReduct (4, d);
    }

    beforeReduct (dim);
    redLLLNTL (korlat->getBasis(), &sqlen, 0.5, dim, prec);
    if (!red->shortestVector())
       std::cout << " shortestVector failed for LLL, delta = 0.5 \n";
    afterReduct (5, d);

    beforeReduct (dim);
    redLLLNTL (korlat->getBasis(), &sqlen, 0.8, dim, prec);
    if (!red->shortestVector())
       std::cout << " shortestVector failed for LLL, delta = 0.8 \n";
    afterReduct (6, d);

    beforeReduct (dim);
    redLLLNTL (korlat->getBasis(), &sqlen, 0.99999, dim, prec);
    if (!red->shortestVector())
       std::cout << " shortestVector failed for LLL, delta = 0.99999 \n";
    afterReduct (7, d);

    beforeReduct (dim);
    redBKZ (korlat->getBasis(), &sqlen, 0.99999, 10, 0, dim, prec);
    if (!red->shortestVector())
       std::cout << " shortestVector failed for BKZ, delta = 0.99999 \n";
    afterReduct (8, d);
}

/**  
*  In this testing loop, we generate `numRep` multipliers `a` and for each one
*/ 
static void testLoop(long numRep) {
    long d;  // dim = dimensions[d].
    Int a;        // The LCG multiplier
    // IntMat basisdual1, basisdual2;
    // basisdual1.SetDims(maxdim, maxdim); // Will be initial triangular dual basis.
    // basisdual2.SetDims(maxdim, maxdim); // m-dual basis passed to LLL, etc.
    for (d = 0; d < numSizes; d++)   // Reset accumulators.
        for (int64_t meth = 0; meth < numMeth; meth++) {
            timer[meth][d] = 0;     sumSq[meth][d] = 0.0;
        }
    totalTime = clock();
    std::cout << "Before the for loop \n";
    for (int64_t r = 0; r < numRep; r++) {
        a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
        std::cout << "Before seta \n";
        korlat->seta(a);
        std::cout << "After seta \n";
        for (d = 0; d < numSizes; d++) {  // Each matrix size
            // korlat->buildDualBasis(dim);
            std::cout << "a = " << a << ",  dim = " << dimensions[d] << "\n";
            // copy(korlat->getDualBasis(), dualbasis1, dim, dim); // Triangular basis.
            reduceBasis(d);
        }
    }
}

static void printResults() {
    long d;
    std::cout << " dim:    ";
    for (d = 0; d < numSizes; d++)
        std::cout << std::setw(8) << dimensions[d] << "  ";
    std::cout << "\n\n";
    for (int meth = 0; meth < numMeth; meth++) {
        std::cout << names[meth] << " ";
        for (d = 0; d < numSizes; d++)
            std::cout << std::setw(9) << timer[meth][d] << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "Sums of square lengths of shortest basis vector:\n";
    std::cout << " dim:    ";
    for (d = 0; d < numSizes; d++)
        std::cout << std::setw(13) << dimensions[d] << "  ";
    std::cout << "\n\n";
    for (int meth = 0; meth < numMeth-2; meth++) {
        std::cout << names[meth] << "  ";
        for (d = 0; d < numSizes; d++)
            std::cout << std::setw(14) << sumSq[meth][d] << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "Total time: "
            << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
            << " seconds\n\n\n";
}

int main() {
    long numRep = 1;    // Number of replications (multipliers) for each case.
    sqlen.SetLength(maxdim);   // Done here because cannot be done in preamble.
    korlat = new Rank1Lattice<Int, Real>(m, maxdim, false, true);
    red = new ReducerBB<Int, Real>(*korlat);

    std::cout << "Types: " << strFlexTypes << "\n";
    std::cout << "PrecisionType: " << prec << "\n\n";
    std::cout << "Results of TestReducersSpeed.cc with m = " << m << "\n";
    std::cout << "Timings (in microseconds) for different methods for "
              << numRep << " replications \n\n";
    testLoop(numRep);
    printResults();
}

/*

//#define TYPES_CODE  LD     // Int == int64_t
#define TYPES_CODE  ZD     // Int == ZZ

#include <iostream>
#include <ctime>  
//#include <NTL/mat_GF2.h>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
//#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
// #include "latticetester/FigureOfMeritM.h"

using namespace LatticeTester;

Int m, a;     // Modulus, LCG Multiplier
// Int m(1048573);  // Prime modulus near 2^{20}
// Int m(1073741827);  // Prime modulus near 2^{30}
// Int m(1099511627791);  // Prime modulus near 2^{40}
// Int m(1125899906842597);  // Prime modulus near 2^{50}

const long numMeth = 6;
const long numSizes = 4;    // Number of matrix sizes (choices of dimension).
const long numPrimes = 2;    // Number of matrix sizes (choices of dimension).
const long dimensions[numSizes] = { 10, 20, 30, 40 };
const long primes[numPrimes] = { 1048573, 1073741827 };
const DecompTypeBB decomp = CHOLESKY; //Decomposition type inside BB algorithm
const bool choiceDual = false; // Choose whether to use primal or dual basis
const long width = 10;

IntMat basis1;
Rank1Lattice<Int, double> *korlat;  // The original Korobov lattice.
IntLattice<Int, double> *lat;       // A copy of the lattice whose basis can be modified.
Reducer<Int, Real> *red;            // The Reducer object that we use.

long dim; //Saves the current dimension

clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes];
clock_t total_times[numSizes];

std::string names[numMeth] = { "None     ", "pairwise ", "LLL5     ", "LLL8     ",
      "LLL99999 ", "BKZ      " };

void copyLattice (const bool & du) {
   // du decides if the primal (false) or dual (true) lattice is used
   if (!du)
      copy(korlat->getBasis(), basis1);
   else
      copy(korlat->getDualBasis(), basis1);
   lat = new IntLattice<Int, Real>(basis1, m, dim);
   red->setIntLattice(*lat);
   lat->updateVecNorm();     
}

void shortestNone (IntLattice<Int, Real> & lattice) {
    // Finding shortest vector without pre-reduction takes too much in big dimensions
   if (dim < 30) {
      if (!red->shortestVector(lattice)) {
         std::cout << " shortestVector failed with no pre-reduction, dim  = " << dim << "\n";
      }
   }    
}

void shortestDieter (IntLattice<Int, Real> & lattice) {
   red->redDieter(0);
   if (dim < 30) {
      if (!red->shortestVector(lattice)) {
         std::cout << " shortestVector failed for pairwise pre-reduction with dim  = " << dim << "\n";
      }
   }
}

void shortestBKZ (IntLattice<Int, Real> & lattice) {
   red->redBKZ(lattice.getBasis());
   if (!red->shortestVector(lattice))
      std::cout << " shortestVector failed for BKZ \n";
}

void shortestLLL (IntLattice<Int, Real> & lattice, double delta) {
   red->redLLLNTL(lattice.getBasis(), delta);
   if (!red->shortestVector(lattice))
      std::cout << " shortestVector failed for LLL " << delta << "\n";  
}

void printOutput() {
   // Create output string
  std::cout << "Results of ReducerComparison.cc with m = " << m << "\n";
   std::cout << "Flexible types: " << strFlexTypes << "\n";
   std::cout << "Timings for finding shortest vector with BB, with various pre-reduction methods.\n";
   std::cout << "The times are in basic clock units. \n";
   std::cout << " dim:    ";
   for (int d = 0; d < numSizes; d++)
      std::cout << std::setw(width) << dimensions[d] << " ";
   std::cout << std::endl << std::endl;
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      for (int d = 0; d < numSizes; d++)
         std::cout << std::setw(width) << timer[meth][d] << " ";
      std::cout << std::endl;
   }
   std::cout << std::endl;
  std::cout << "Total time: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC) << " seconds\n\n\n";
 
}


int main() {
  lat = new IntLattice<Int, Real>(m, dimensions[numSizes]);
  red = new Reducer<Int, Real>(*lat);
  red->setDecompTypeBB(decomp);
  for (int p = 0; p < numPrimes; p++) {
     totalTime = clock();
     //Go through the list of prime numbers and set a accordingly
     m = primes[p];
     a = m/7;
     for (int d = 0; d < numSizes; d++) {
  
        dim = dimensions[d];
        korlat = new Rank1Lattice<Int, Real>(m, a, dim, choiceDual);
        korlat->buildBasis(dim);
        basis1.SetDims(dim, dim);
        
        copyLattice(choiceDual);   
        tmp = clock();
        shortestNone(*lat);
        timer[0][d] = clock() - tmp;
     
        copyLattice(choiceDual);
        tmp = clock();
        red->redDieter(0);
        shortestDieter(*lat);
        timer[1][d] = clock() - tmp;
        
        copyLattice(choiceDual);
        tmp = clock();
        shortestLLL(*lat, 0.5);
        timer[2][d] = clock() - tmp;
   
        copyLattice(choiceDual);
        tmp = clock();
        shortestLLL(*lat, 0.8);
        timer[3][d] = clock() - tmp;
      
        copyLattice(choiceDual);
        tmp = clock();
        shortestLLL(*lat, 0.9999);
        timer[4][d] = clock() - tmp;
     
        copyLattice(choiceDual);
        tmp = clock();
        shortestBKZ(*lat);
        timer[5][d] = clock() - tmp;
     }       
     printOutput(); 
  }
  return 0;
}


Results of ReducerComparison.cc with m = 1048573
Flexible types: Int = NTL::ZZ, Real = double
Timings for finding shortest vector with BB, with various pre-reduction methods.
The times are in basic clock units. 
 dim:            10         20         30         40 

None           71331     572928          1          1 
pairwise        3010     237591       7650      17777 
LLL5             182       2003     246096     701089 
LLL8             164       1408      55743      97849 
LLL99999         176       1522      26298      67121 
BKZ              318       2566      23969      48528 

Total time: 2.18838 seconds

Results of ReducerComparison.cc with m = 1073741827
Flexible types: Int = NTL::ZZ, Real = double
Timings for finding shortest vector with BB, with various pre-reduction methods.
The times are in basic clock units. 
 dim:            10         20         30         40 

None         2649188    6738802          1          1 
pairwise      112822    2262319       8251      18788 
LLL5             219       9158     545774  250225757 
LLL8             239       1893     104425   17262999 
LLL99999         263       2383      19737    3125776 
BKZ              458       3894      23967    3872037 

Total time: 286.992 seconds

*/