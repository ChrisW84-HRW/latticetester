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
#include "latticetester/IntLattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"


namespace LatticeTester {

/**
 * The test is based on the known shortest vector lengths of some 
 * low dimensional lattices as found in F. Pausinger: 'Long shortest vectors 
 * in low dimensional lattices', Discrete Mathematics, 12 (2020).
 * The BB algorithm is tested against these vector lenghts and 
 * is supposed to yield the correct results. 
**/
template<typename Int, typename Real>
void testRedBB() {
    
    // Test Case 1 (according to Section 3)
    IntMat base;
    Int N(288);
    base.SetDims(3, 3);
    base[0][0] = 19; base[0][1] = 3; base[0][2] = 39;
    base[1][0] = 3; base[1][1] = 39; base[1][2] = 19;
    base[2][0] = -39; base[2][1] = -19; base[2][2] = 3;
    IntLattice<Int, Real> lat(base, N, 3);
    Real m_minBB;
    ReducerBB<Int, Real> red(4);  
    if (!red.shortestVector(lat)) m_minBB = 0;
    else m_minBB = red.getMinLength2(); 
    CHECK(m_minBB == 1891);

    // Test Case 2 (following from Theorem 2)
    Int m(4103);
    N = power(m,3)-1;
    base[0][0] = 0; base[0][1] = 0; base[0][2] = N;
    base[1][0] = 0; base[1][1] = N; base[1][2] = 0;
    base[2][0] = 1; base[2][1] = m; base[2][2] = power(m,2);
    IntLattice<Int, Real> lat2(base, N, 3);
    if (!red.shortestVector(lat2)) m_minBB = 0;
    else m_minBB = red.getMinLength2(); 
    CHECK(m_minBB == NTL::to_RR(1 + power(m,2) + power(m,4)));

    // Test Case 3 (following from Theorem 4)
    m = 371;
    Int k;
    k = (m-1)/2;
    base.SetDims(4, 4);
    N = (power(m, 4) - 1)/2;
    base[0][0] = 0; base[0][1] = 0; base[0][2] = 0; base[0][3] = N;
    base[1][0] = 0; base[1][1] = 0; base[1][2] = N; base[1][3] = 0;
    base[2][0] = 0; base[2][1] = N; base[2][2] = 0; base[2][3] = 0;
    base[3][0] = 1; base[3][1] = m; base[3][2] = power(m,2); base[3][3] = power(m,3);
    IntLattice<Int, Real> lat3(base, N, 4);
    if (!red.shortestVector(lat3)) m_minBB = 0;
    else m_minBB = red.getMinLength2(); 
    CHECK(m_minBB == NTL::to_RR(4 * k * k * (4*k*k+6*k+3)*(2*k*(k+1)+1)));       
}


TEST_CASE("Test reducer static") {
     
     testRedBB<NTL::ZZ, double> ();
     testRedBB<NTL::ZZ, NTL::RR> ();
}

};
