/**
 * This example illustrates the usage of the BasisConstruction module.
 *
 * TO REWRITE WHEN DONE:   *****
 *
 * GOALS:
 *
 * -- perform speed comparisons between the methods to construct a basis or m-dual basis
 *    for a general lattice, given a set of generating vectors for the primal.
 *    We want to do this for various numbers of dimensions, as done below.
 * -- comparing the speed with ZZ vs int64, with double vs RR, etc.
 * -- Show how much overhead there is when we resize the IntMat objects frequently.
 *
 *  ================
 * This reads matrices from files and builds a basis and a dual for an `IntLattice`
 * object. The files this is set to use are in the `bench.zip` archive. To
 * execute the program, the archive should be unziped and the `bench` folder
 * should be put in the same directory from which the executable is called.

 * The bases we used to showcase of BasisConstruction methods are in a folder named 
 * 'examples/bench/'. Each file in 'examples/bench/' folder contain a basis, and the 
 * file is nameed as follows: 'prime_dimBasis_exanpleNumber' where 'prime' is modulo 
 * value of the basis, 'dimBasis' is the dimension of the basis, and 'exampleNumber' 
 * is the number of the example for the bases of dimension 'dimBasis'.
 *
 * This example reads matrices from files and performs the different construction
 * algorithms in BasisConstruction on them. The program then prints the execution
 * time of the various algorithms. Note that the execution of the program is not
 * what you would expect in reality since bench contains random full matrices.
 *
 * We show a use of BasisContruction::upperTriangularBasis, 
 * BasisContruction::lowerTriangularBasis, BasisContruction::LLLConstruction with 
 * two differrent parametter of delta, BasisContruction::mDualUpperTriangular, 
 *  BasisContruction::mDualBasis
 *
 * In this example, we can compare the speed of BasisConstruction<Int>::calcDual method
 * which compute an m-dual basis using any basis in input,
 * and BasisConstruction::mDualUpperTriangular method which compute an m-dual basis
 * with an upper triangular basis.
 *
 * We can also compare the speed of 'BasisConstruction::upperTriangularBasis'
 * and the speed of 'BasisConstruction::LLLConstruction'
 *
 * Example of results with m = 1048573 (prime modulus near 2^{20}):
 *
 *  ***  These must be replaced  ***
 *
 **/

// Select the flexible types Int and Real here.
//#define TYPES_CODE  LD     // Int == int64_t
//#define TYPES_CODE  ZD     // Int == ZZ, Real = double
#define TYPES_CODE  ZR     // Int == ZZ, Real = RR

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
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/ReducerStatic.h"
// #include "latticetester/LLL_FPInt.h"

using namespace LatticeTester;

//Int m(1021);     // Modulus m = 1021
Int m(1048573);  // Prime modulus near 2^{20}
//Int m(1073741827);  // Prime modulus near 2^{30}
//Int m(1099511627791);  // Prime modulus near 2^{40}
//Int m(1125899906842597);  // Prime modulus near 2^{50}
Int a;       // The LCG multiplier

const long numSizes = 5;    // Number of matrix sizes (choices of dimension).
//const long dimensions[numSizes] = { 10, 12, 14, 16, 18, 20 };
const long dimensions[numSizes] = { 4, 6, 10, 20, 30 };
long maxdim = dimensions[numSizes - 1];   // Maximum dimension
const long numMeth = 6;    // Number of methods, and their names.
std::string names[numMeth] = { "LLL5      ", "LLL9      ", "LLL99999  ",
        "LLL99999-R", "UppTri    ", "mDualUT   "};
// PrecisionType prec = QUADRUPLE;  // QUADRUPLE, XDOUBLE, RR
// PrecisionType prec = DOUBLE;

// Here we use ctime directly for the timings, to minimize overhead.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes];
Real sumSq[numMeth][numSizes];
RealVec sqlen; // The length of this vector is set in the main.

// Run speed test for dim = dimensions[d], with given matrices.
static void transformBases (long d, long dim, IntMat &basis1, IntMat &basis2,
        IntMat &basisdual) {
    // We apply LLL to basis1 with different values of `delta`, incrementally.
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    tmp = clock();
    LLLConstruction0(basis2, 0.5, dim, dim, &sqlen);
    timer[0][d] += clock() - tmp;
    sumSq[0][d] += sqlen[0];

    // We continue the LLL process with a larger `delta`.
    tmp = clock();
    LLLConstruction0(basis2,  0.9, dim, dim, &sqlen);
    timer[1][d] += clock() - tmp;
    sumSq[1][d] += (sqlen)[0];

    tmp = clock();
    LLLConstruction0(basis2, 0.99999, dim, dim, &sqlen);
    timer[2][d] += clock() - tmp;
    sumSq[2][d] += (sqlen)[0];

    // Here we restart LLL from the initial triangular basis.
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    tmp = clock();
    LLLConstruction0(basis2, 0.99999, dim, dim, &sqlen);
    timer[3][d] += clock() - tmp;
    sumSq[3][d] += (sqlen)[0];

    // We now construct an upper-triangular basis from basis2 into basis1.
    tmp = clock();
    upperTriangularBasis(basis2, basis1, m, dim, dim);
    timer[4][d] += clock() - tmp;

    // We compute an m-dual basis to basis1.
    tmp = clock();
    mDualUpperTriangular(basis1, basisdual, m, dim);
    timer[5][d] += clock() - tmp;
}

/*
// A verbose version, for checking.
static void transformBasesVerb (long d, long dim, IntMat &basis1, IntMat &basis2,
        IntMat &basisdual) {
    // We apply LLL to basis1 with different values of `delta`, incrementally.
    Int sqlength;
    std::cout << "\n ******** dim = " << dim << "\n";

    // ****   Use CopyPartMat instead!
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    //  copy(basis1, basis2, dim, dim);
    std::cout << "Before LLL 0.5: sqlen = " << sqlen << "\n";
    tmp = clock();
    // redLLLNTL(basis2, 0.5, dim, &sqlen);
    LLLConstruction0(basis2, &sqlen, 0.5, dim, dim);   // FAILS for RR
    timer[0][d] += clock() - tmp;
    sumSq[0][d] += sqlen[0];
    std::cout << "After LLL 0.5: sqlen[0] = " << (sqlen)[0] << "\n";
    std::cout << "First basis vector: " << basis2[0] << "\n";
    ProdScal<Int>(basis2[0], basis2[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // We continue the LLL process with a larger `delta`.
    // copy(basis1, basis2, dim, dim);
    tmp = clock();
    // redLLLNTL(basis2, 0.9, dim, &sqlen);
    LLLConstruction0(basis2, &sqlen, 0.9, dim, dim);
    timer[1][d] += clock() - tmp;
    //sumSq[1][d] += sqlen[0];
    std::cout << "After LLL 0.9:  sqlen[0] = " << (sqlen)[0] << "\n";
    ProdScal<Int>(basis2[0], basis2[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // copy(basis1, basis2, dim, dim);
    tmp = clock();
    // redLLLNTL(basis2, 0.99999, dim, &sqlen);
    timer[2][d] += clock() - tmp;
    //sumSq[2][d] += sqlen[0];
    std::cout << "After LLL 0.99999:  sqlen[0] = " << (sqlen)[0] << "\n";
    std::cout << "First basis vector: " << basis2[0] << "\n";
    ProdScal<Int>(basis2[0], basis2[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // Here we restart LLL from the initial triangular basis.
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    // copy(basis1, basis2, dim, dim);
    tmp = clock();
    LLLConstruction0(basis2, &sqlen, 0.99999, dim, dim);
    timer[3][d] += clock() - tmp;
    //sumSq[3][d] += sqlen[0];

    // We now construct an upper-triangular basis from basis2 into basis1.
    tmp = clock();
    upperTriangularBasis(basis2, basis1, m, dim, dim);
    timer[4][d] += clock() - tmp;

    // We compute an m-dual basis to basis1.
    tmp = clock();
    mDualUpperTriangular(basis1, basisdual, m, dim);
    timer[5][d] += clock() - tmp;

    // mDualBasis is currently implemented only for Int = ZZ and dim = maxDim.
    // BasisConstruction<Int>::mDualBasis(basis2, basisdual, m);
}
*/

// In this testing loop, new `Rank1Lattice` objects are created
// and the  `IntMat` matrices are resized inside the loop.
static void testLoopResize(long numRep) {
    long d, dim;
    IntMat basis1, basis2, basisdual;
    Rank1Lattice<Int, Real> *korlat;    // Will be a Korobov lattice.
    for (d = 0; d < numSizes; d++)      // Reset timers and sums.
        for (int64_t meth = 0; meth < numMeth; meth++) {
            timer[meth][d] = 0;   sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
        }
    totalTime = clock();
    for (int64_t r = 0; r < numRep; r++) {
        a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
        //for (d = 0; d < 2; d++) {  // Each matrix size
        for (d = 0; d < numSizes; d++) {  // Each matrix size
            dim = dimensions[d]; // The corresponding dimension.
            basis1.SetDims(dim, dim); // Will be initial triangular basis.
            basis2.SetDims(dim, dim); // Will be LLL-reduced basis.
            basisdual.SetDims(dim, dim);  // m-dual basis.
            korlat = new Rank1Lattice<Int, Real>(m, a, dim, true, true);
            korlat->buildBasis(dim);
            basis1 = korlat->getBasis();
            transformBases(d, dim, basis1, basis2, basisdual);
            delete korlat;
        }
    }
    basis1.kill();  // Since we create objects repeatedly,
    basis2.kill();  // it is a good idea to release the memory when we are done.
    basisdual.kill();
}

// In this testing loop, we try to minimize the creation of objects.
// The `IntMat` and `Rank1Lattice` objects are created only once.
static void testLoopNoResize(long numRep) {
    long d, dim;  // Index of dimension.
    IntMat basis1, basis2, basisdual;
    basis1.SetDims(maxdim, maxdim); // Will be initial triangular basis.
    basis2.SetDims(maxdim, maxdim); // Will be LLL-reduced basis.
    basisdual.SetDims(maxdim, maxdim);  // m-dual basis.
    Rank1Lattice<Int, Real> *korlat;  // We create a single Korobov lattice object.
    korlat = new Rank1Lattice<Int, Real>(m, maxdim, true, true);

    for (d = 0; d < numSizes; d++)   // Reset accumulators.
        for (int64_t meth = 0; meth < numMeth; meth++) {
            timer[meth][d] = 0;     sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
        }
    totalTime = clock();
    for (int64_t r = 0; r < numRep; r++) {
        a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
        korlat->seta(a);
        for (d = 0; d < numSizes; d++) {  // Each matrix size
            dim = dimensions[d]; // The corresponding dimension.
            korlat->buildBasis(dim);
            std::cout << "a = " << a << ",  dim = " << dimensions[d] << "\n";
            copy(korlat->getBasis(), basis1, dim, dim); // Triangular basis.
            transformBases(d, dim, basis1, basis2, basisdual);
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
    std::cout << "Sums of square lengths of shortest basis vector";
    std::cout << " (must be the same across all implementations):\n";
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
    long numRep = 1000;    // Number of replications (multipliers) for each case.
    sqlen.SetLength(maxdim);   // Done here because cannot be done in preamble.
    std::cout << "Types: " << strFlexTypes << "\n";
    std::cout << "PrecisionType: " << prec << "\n\n";
    std::cout << "Results of TestBasisConstructionSpeed.cc with m = " << m << "\n";
    std::cout << "Timings for different methods, in basic clock units \n\n";
    testLoopNoResize(numRep);
    std::cout << "Timings for `testLoop No Resize`\n";
    printResults();
    testLoopResize(numRep);
    std::cout << "Timings for `testLoop Resize` (many objects are created or resized)\n";
    printResults();

}

