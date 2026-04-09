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
#include "latticetester/Weights.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"

namespace LatticeTester {

/**
 * Takes a modulus 'm' and a multiplier 'a' and checks if the
 * figure of merit for the primal and dual lattice does not depend
 * on what type of pre-reduction is used. Futhermore it checks that
 * the result does not change depending if we call the FoMs for the
 * succssive and non-succesive coordinates seperately and take their
 * minimum.
 */
template<typename Int, typename Real>
void fomTest(Int m, Int a, int64_t maxdim, int dim, NTL::Vec < int64_t > t) {
     Rank1Lattice<Int, Real> lat(m, maxdim);
     IntLattice<Int, Real> proj(m, maxdim);
     string s;
     NormType norm = L2NORM;
     WeightsUniform weights(1.0);
     lat.seta(a);
     lat.buildBasis(dim);
     lat.buildDualBasis(dim);
     NormaBestLat normaPrimal(log(m), 1, maxdim, norm);  // Factors computed for primal. 
     NormaBestLat normaDual(-log(m), 1, maxdim, norm);  // Factors computed for dual. 
     ReducerBB<Int, Real> red(maxdim);   // Single ReducerBB with internal lattice `lat`.
     FigureOfMeritM<Int, Real> fomPrimal(t, weights, normaPrimal, &red, true); // includeFirst = true.
     FigureOfMeritDualM<Int, Real> fomDual(t, weights, normaDual, &red, true); // FoM for dual lattice.
     double merit1, merit2, merit3, merit4, minmerit;
     
     // Tests for the primal lattice
     
     // Use ususal method with BKZ prereduction
     merit1 = fomPrimal.computeMerit(lat, proj);
     // Use LLL instead of BKZ as prereduction
     fomPrimal.setBKZ(0);
     fomPrimal.setLLL(0.99);
     merit2 = fomPrimal.computeMerit(lat, proj);
     CHECK(merit1 == merit2);
     // No prereduction
     fomPrimal.setBKZ(0);
     fomPrimal.setLLL(0);
     merit2 = fomPrimal.computeMerit(lat, proj);
     CHECK(merit1 == merit2);

     // Now calculate FoM for successive and non-successive coordinates separately
     merit3 = fomPrimal.computeMeritNonSucc(lat, proj);
     merit4 = fomPrimal.computeMeritSucc(lat);
     minmerit = min(merit3, merit4);
     CHECK(merit1 == minmerit);

     // Tests for the dual lattice

     merit1 = fomDual.computeMerit(lat, proj);
     // Use LLL instead of BKZ as prereduction
     fomDual.setBKZ(0);
     fomDual.setLLL(0.99);
     merit2 = fomDual.computeMerit(lat, proj);
     CHECK(merit1 == merit2);
     // No prereduction
     fomDual.setBKZ(0);
     fomDual.setLLL(0);
     merit2 = fomDual.computeMerit(lat, proj);
     CHECK(merit1 == merit2);
   
     // Now calculate FoM for successive and non-successive coordinates separately
     merit3 = fomDual.computeMeritNonSucc(lat, proj);
     merit4 = fomDual.computeMeritSucc(lat);
     minmerit = min(merit3, merit4);
     CHECK(merit1 == minmerit);   
       
}


TEST_CASE("Test figure of merit") {
     int64_t maxdim = 8;
     int dim = 8;
     NTL::ZZ m; m = 1021;
     NTL::ZZ a; a = 822;
     NTL::Vec < int64_t > t;
     t.SetLength(4);
     t[0] = 4;    
     t[1] = 8;    // For pairs.
     t[2] = 6;    // For triples.
     t[3] = 4;    // For quadruples.

     fomTest<NTL::ZZ, double> (m, a, maxdim, dim, t);
     fomTest<NTL::ZZ, NTL::RR> (m, a, maxdim, dim, t);

     t.SetLength(1);
     t[0] = 2;
     
     a = 703;
     fomTest<NTL::ZZ, double> (m, a, maxdim, dim, t);
     fomTest<NTL::ZZ, NTL::RR> (m, a, maxdim, dim, t);
}

};
