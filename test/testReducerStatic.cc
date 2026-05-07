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
 * Takes a modulus 'm' and a multiplier 'a' and checks if the
 * if the shortes vectors after LLL and BKZ reduction satisfy 
 * the known theoretical inequalities as reported in the gudie.
 * The values for delta are delta = 0.5 and delta = 0.999. 
 * The blocksize in BKZ is equal to 'k'.
 */
template<typename Int, typename Real>
void testRedStat(Int m, Int a, int64_t maxdim, int64_t k) {
     Rank1Lattice<Int, Real> lat(m, maxdim);
     IntLattice<Int, Real> proj(m, maxdim);
     lat.seta(a);
     Real m_minLLL, m_minBKZ, m_minBB;
     ReducerBB<Int, Real> red(maxdim);  

     // Test if the shortest vector length after LLL satisfies the theoretically known inequality for delta = 0.5
     double delta(0.5);
     lat.buildBasis(maxdim);
     m_minLLL = redLLL<Int, Real>(lat.getBasis(), delta, 0);
     if (!red.shortestVector(lat)) m_minBB = 0;
     else m_minBB = red.getMinLength2(); 
     CHECK(m_minLLL/m_minBB < pow(delta-1/4,(1-maxdim)/2));

     // Do the same for delta = 0.999
     delta = 0.999;
     lat.buildBasis(maxdim);
     m_minLLL = redLLL<Int, Real>(lat.getBasis(), delta, 0);
     if (!red.shortestVector(lat)) m_minBB = 0;
     else m_minBB = red.getMinLength2(); 
     CHECK(m_minLLL/m_minBB < pow(delta-1/4,(1-maxdim)/2));

     // Now check the inequality for the BKZ algorithm and delta = 0.5     
     double ak;
     ak = pow(k,1+log(k));
     delta = 0.5;
     lat.buildBasis(maxdim);
     m_minBKZ = redBKZ<Int, Real>(lat.getBasis(), delta, k);
     if (!red.shortestVector(lat)) m_minBB = 0;
     else m_minBB = red.getMinLength2(); 
     CHECK(m_minBKZ/m_minBB < pow(ak,(maxdim-1)/(2*(k-1))));

     // and for delta = 0.999;
     delta = 0.999;
     lat.buildBasis(maxdim);
     m_minBKZ = redBKZ<Int, Real>(lat.getBasis(), delta, k);
     if (!red.shortestVector(lat)) m_minBB = 0;
     else m_minBB = red.getMinLength2(); 
     CHECK(m_minBKZ/m_minBB < pow(ak,(maxdim-1)/(2*(k-1))));
       
}


TEST_CASE("Test reducer static") {
     int64_t maxdim = 31;
     NTL::ZZ m(1048573);
     NTL::ZZ a(82281);
     int64_t k;
     k = (maxdim-1)/2;
     
     testRedStat<NTL::ZZ, double> (m, a, maxdim, k);
     testRedStat<NTL::ZZ, NTL::RR> (m, a, maxdim, k);
}

};
