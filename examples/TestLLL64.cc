/**
 * This example illustrates the use of functions from the `BasisConstruction` class,
 * with a small five-dimensional lattice obtained from an LCG with a small modulus.
 * The basis is printed between function calls, to show what is going on.
 * An initial upper-triangular basis is constructed by `Rank1Lattice` in a standard way.
 * We then call LLL with `delta = 0.5` to obtain a basis with shorter vectors,
 * which is not triangular. Calling `upperTriangulatBasis` transforms this basis
 * to an upper triangular basis, which happens to be the same as the initial one.
 * Different triangularization methods are compared.  Then LLL with different
 * values of `delta`. Then a dual basis is computed in different ways.
 *
 * After that, we look at projections of this lattice over subsets of coordinates.
 * We show how to construct a basis for such a projection, compute the corresponding
 * dual basis, and compute a shortest vector in this dual basis.
 *
 * This program can be run with the two different types of `Int` (just change the first
 * line below), and with various choices of the modulus `m` and multiplier `a`.
 **/

// define TYPES_CODE  LD     // int64_t
#define TYPES_CODE  ZD     // ZZ

#include <iostream>
#include <cstdint>
// #include <ctime>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
//#include "latticetester/Reducer.h"
#include "latticetester/Coordinates.h"

using namespace LatticeTester;

//Int m(101);      // Modulus m = 101
//Int m(1021);     // Modulus m = 1021
Int m(1048573);  // Modulus m = 1048573 (prime number near 2^{20})
Int a(209714);       // The LCG multiplier
const long dim = 6;  // Dimension of lattice.
const long maxDim = 30;  // Dimension of lattice.

int main() {
    std::cout << "Types: " << strFlexTypes << "\n\n";
    std::cout << "TestLLL64 \n";

    // All the IntMat objects are created in 10 dimensions, but we may use less.
    IntMat basis1, basis2;
    basis1.SetDims(maxDim, maxDim);
    basis2.SetDims(maxDim, maxDim);
    long sqlength;

    // We construct a Korobov lattice in dim dimensions.
    Rank1Lattice<Int, double> *korlat;
    korlat = new Rank1Lattice<Int, Real>(m, a, maxDim, true, true);
    korlat->buildBasis(dim);
    copy(korlat->getBasis(), basis1, dim, dim);  // This initial basis is triangular.
    std::cout << "Initial Korobov lattice basis = \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    BasisConstruction<Int>::LLLConstruction0(basis1, 0.99999);
    std::cout << "Basis after LLL with delta=0.99999: \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";
    return 0;
}
