#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <iostream>
#include <cstdint>
#include <cmath>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Coordinates.h"

namespace LatticeTester {

/**
 * Tests different methods of BasisConstruction.h for 
 * the modulus 'm' and the dimension 'maxdim'. It checks that
 * the functions 'lowerTriangularBasis" and 'upperTriangularBasis'
 * indeed yield lower and upper triangular basis matrices.
 * Moreover, it is checked that the three methods to calculate the 
 * m-dual matrix indeed yield the dual of the current base.
 * Finally, the construction of projections is checked.
**/
template<typename Int, typename Real>
void BasisConstructionTest(Int m, int64_t maxdim) {
     IntMat basis, basis2, basis3;
     basis.SetDims(maxdim, maxdim);
     basis2.SetDims(maxdim, maxdim);
     basis3.SetDims(maxdim, maxdim);
     // Define a full rank matrix do be used as basis matrix
     for (int i = 0; i < maxdim; i++)
     { 
          for (int j = 0; j < maxdim; j++)
          {
               basis[i][j] = (i+1)*(j+1) + (j-i) + (i==j);
          }
     }
     // Make sure that the determininant is divisible by m.
     for (int j = 0; j < maxdim; j++)
         basis[0][j] *= m;
     basis2 = basis;
          
     IntLattice<Int, Real> lat(basis, m, maxdim);     
     // Check that lower and upper triangular basis construction really produce lower / upper triangular basis matrices
     // At the same time check the corresponding functions to build the dual
     lowerTriangularBasis<Int> (lat.getBasis(), basis, m);
     for (int i = 0; i < maxdim; i++)
     { 
          for (int j = i + 1; j < maxdim; j++)
          {
               CHECK(lat.getBasis()[i][j] == 0);
          }
     }
     mDualLowerTriangular<Int> (lat.getDualBasis(), lat.getBasis(), m);
     CHECK(lat.checkDuality() == true);  
     
     basis = basis2;
     upperTriangularBasis<Int> (lat.getBasis(), basis, m);
     for (int i = 0; i < maxdim; i++)
     { 
          for (int j = 0; j < i; j++)
          {
               CHECK(lat.getBasis()[i][j] == 0);
          }
     }
     mDualUpperTriangular<Int> (lat.getDualBasis(), lat.getBasis(), m);
     CHECK(lat.checkDuality() == true);  

     basis = basis2;
     lat.setBasis(basis, maxdim);
     mDualBasis<Int>(lat.getDualBasis(), basis, m);
     CHECK(lat.checkDuality() == true);  
     
     // Check if projections function correctly
     Coordinates c;
     int bound;
     bound = maxdim / 2;
     // The coordinates are set to contain every second dimension
     for (int j = 0; j < bound; j++)
         c.insert(2*j+1);     
     basis = basis2;
     projectMatrix(basis3, basis, c);     
     for (int i = 0; i < maxdim; i++) {
          for (int j = 0; j < bound; j++)
              CHECK(basis3[i][j] == basis[i][2*j]);
     }
}

TEST_CASE("Test basis construction") {
     int64_t maxdim = 8;
     NTL::ZZ m; m = 1021;
     BasisConstructionTest<NTL::ZZ, double>(m, maxdim);
};

}
