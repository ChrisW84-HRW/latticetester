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

#include "NTL/LLL.h"

#include "latticetester/ntlwrap.h"
#include "latticetester/Util.h"

namespace LatticeTester {
  /**
   * **WARNING** THIS MODULE IS UNDER CONSTRUCTION.
   *
   * This class implements general methods to perform a lattice basis
   * construction from a set of vectors, as well as general methods to obtain
   * the dual lattice basis depending on the current lattice basis.
   *
   * This module still has to be implemented.
   *
   * This module should work only with NTL matrices because those are objects
   * aware of their shape. I think that this is important since we use object
   * oriented stuff.
   *
   * What is done here seems like very simple linear algebra BUT it is not. The
   * fact that we operate only on the ring of integers makes it so that these
   * simple algorithms can become way too slow very rapidly because the numbers
   * they manipulate grow very fast. This also means that in many cases the
   * usage of standard `long` type integers will overflow.
   * */
  template<typename BasInt> class BasisConstruction{

    private:
      typedef NTL::vector<BasInt> BasIntVec;
      typedef NTL::matrix<BasInt> BasIntMat;

    public:
      /**
       * */
      BasisConstruction(){};

      /**
       * Destructor.
       * */
      ~BasisConstruction(){};

      /**
       * This functions takes a set of generating vectors of a vector space and
       * finds a basis for this space whilst applying LLL reduction. This is
       * much faster than applying GCDConstruction, but it doesn't help building
       * the dual.
       * */
      void LLLConstruction(BasIntMat& matrix);

      /**
       * This function does some kind of Gaussian elimination on
       * \f$\text{span}(\text{matrix}^t)\f$ as a vector space on \f$\mathbb{Z}\f$ to
       * obtain a basis of this vector space (lattice).
       * This takes a matrix `matrix` where each row `v[0], ..., v[n]` is
       * considered as a vector. This basically performs, for each column,
       * Euclid's algorithm on the rows under the diagonal to set all
       * coefficients to zero except the one in the diagonal. Finding the GCD
       * of all the coefficients from the diagonal and under allows to perform
       * Gaussian elimination using only the 3 allowed operations that do not
       * change the span of the vectors:
       *   - Multiply a line by \f$-1\f$,
       *   - Add an integer multiple of line \f$i\f$ to line \f$j\f$ for \f$i \neq j\f$,
       *   - Swap line \f$i\f$ and line \f$j\f$.
       * After constructing this basis, the algorithm eliminates negative
       * coefficients in the matrix.
       *
       * WATCH OUT. This function (and building mecanism) are very memory heavy
       * has the numbers below the diagonal can grow very big.
       *
       * \todo Give a running time. **Marc-Antoine**: This runs in
       * O(min(row, col)*lg(n)*row + row^2) time I think, but it might be a bit
       * more because we do operations on vectors so we might need to multiply
       * this by col.
       * */
      void GCDConstruction(BasIntMat& matrix);

      /**
       * Suppose the basis matrix contains basis vectors on its lines and is
       * \f$p\times q\f$ where \f$q \geq p\f$. We can compute the \f$m\f$-dual
       * as follows.
       * Let's note the basis matrix `V` and the dual matrix `W` and have lines
       * of `W` also contain dual basis vectors in its lines. We know that
       * \f$VW^t = mI_{p\times p}\f$ where \f$I\f$ is the identity matrix. Now, in
       * the case of a \f$m\f$-dual computation, we can assume that all the
       * arithmetic is done modulo \f$m\f$ and that \f$m\f$ is prime??
       * */
      void DualConstruction(BasIntMat& matrix, BasIntMat& dualMatrix,
          BasInt& modulo);

      void DualConstruction2(BasIntMat& matrix, BasIntMat& dualMatrix,
          BasInt& modulo);

    private:

      /**
       * This is the matrix this object is working on. This matrix is modified
       * as the program advances.
       * */
      BasIntMat m_mat;

      /**
       * This stores the matrix the object was initialized with. This is usefull
       * for implementing methods that verify this class does not change the
       * lattice. This is also usefull to compare different reduction methods.
       * Finally, this can be kept for the sake of comparison.
       * */
      BasIntMat m_old;

  };

  //============================================================================
  // Implementation

  template<typename BasInt>
    void BasisConstruction<BasInt>::LLLConstruction(BasIntMat& matrix){
      long rank = NTL::LLL_XD(matrix);
      long num = matrix.NumRows();
      for (long i = 0; i < rank; i++) {
        NTL::swap(matrix[i], matrix[num-rank+i]);
      }
      matrix.SetDims(rank, matrix.NumCols());
    }

  template<typename BasInt>
    void BasisConstruction<BasInt>::GCDConstruction(BasIntMat& matrix)
  {
    // It is important to note that the lines of matrix are the basis vectors
    long rows = matrix.NumRows();
    long cols = matrix.NumCols();
    long max_rank = rows < cols ? rows : cols;
    long rank = 0;
    // The basis will have at most max_rank vectors.
    BasInt q;
    for (long i = 0; i < max_rank; i++) {
      // We find gcd(matrix[i][i], ..., matrix[rows-1][i]) using Euclid
      // algorithm and applying transformations to the matrix
      for (long j = i+1; j < rows; j++) {
        while (matrix[j][i] != 0) {
          matrix[i].swap(matrix[j]);
          q = matrix[j][i]/matrix[i][i];
          matrix[j] -= q * matrix[i];
          //std::cout << "Matrix:\n" << matrix;
        }
      }
      // We make sure that the coefficients are positive.
      // This is because the algorithms we use work for positive vectors.
      // if (matrix[i][i] < 0) matrix[i] *= -1;
      // for (long j = i-1; j >= 0; j--) {
      //   if (matrix[j][i] < 0) {
      //     if (-matrix[j][i] > matrix[i][i]){
      //       q = -matrix[j][i]/matrix[i][i] + 1;
      //       matrix[j] += q * matrix[i];
      //     } else {
      //       matrix[j] += matrix[i];
      //     }
      //   }
      // }
      if (matrix[i][i] != 0) rank++;
    }
    //std::cout << "Matrix before truncating\n" << matrix;
    // We remove zero vectors from the basis.
    matrix.SetDims(rank, cols);
    //std::cout << matrix << std::endl;
  }

  //============================================================================

  /**
   * This algorithm calculates the dual as well as the `m` used for rescalling.
   * It checks if this `m` divides the modulo given to the algorithm.
   * Right now, this assumes the basis is triangular, might need to change it.
   * */
  template<typename BasInt>
    void BasisConstruction<BasInt>::DualConstruction(BasIntMat& matrix,
        BasIntMat& dualMatrix, BasInt& modulo)
  {
    // We need to have a triangular basis matrix
    if (! CheckTriangular(matrix, matrix.NumRows(), modulo))
      GCDConstruction(matrix);
    long dim = matrix.NumRows();
    if (dim != matrix.NumCols()) {
      std::cout << "matrix has to be square, but dimensions do not fit.\n";
      return;
    }
    BasInt m;
    m = 1;
    NTL::ident(dualMatrix, dim);
    BasInt gcd;
    for (long i = dim-1; i>=0; i--) {
      for (long j = i+1; j < dim; j++) {
        dualMatrix[i] -= matrix[i][j] * dualMatrix[j];
        matrix[i] -= matrix[i][j] * matrix[j];
      }
      gcd = matrix[i][i];
      for (long j = i; j < dim; j++) {
        gcd = NTL::GCD(gcd, dualMatrix[i][j]);
      }
      gcd = matrix[i][i] / gcd;
      if (gcd != 1) {
        dualMatrix *= gcd;
        m *= gcd;
      }
      for (long j = i; j < dim; j++) {
        dualMatrix[i][j] /= matrix[i][i];
      }
      matrix[i][i] = 1;
    }
  }

  template<typename BasInt>
    void BasisConstruction<BasInt>::DualConstruction2(BasIntMat& matrix,
        BasIntMat& dualMatrix, BasInt& modulo)
  {
    // We need to have a triangular basis matrix
    if (! CheckTriangular(matrix, matrix.NumRows(), BasInt(0)))
      GCDConstruction(matrix);
    long dim = matrix.NumRows();
    if (dim != matrix.NumCols()) {
      std::cout << "matrix has to be square, but dimensions do not fit.\n";
      return;
    }
    if (modulo < 1) {
      std::cerr << "modulo has to be a positive integer.\n";
      return;
    }
    dualMatrix.SetDims(dim, dim);
    for (int i = 0; i < dim; i++) {
      for (int j = i + 1; j < dim; j++)
        NTL::clear (dualMatrix(i,j));

      if (!NTL::IsZero(matrix(i,i))) {
        BasInt gcd = NTL::GCD(modulo, matrix(i,i));
        modulo *= matrix(i,i) / gcd;
        dualMatrix *= matrix(i,i) / gcd;
      }

      DivideRound (modulo, matrix(i,i), dualMatrix(i,i));

      for (int j = i - 1; j >= 0; j--) {

        NTL::clear (dualMatrix(i,j));

        for (int k = j + 1; k <= i; k++)
          dualMatrix(i,j) += matrix(j,k) * dualMatrix(i,k);

        if (dualMatrix(i,j) != 0)
          dualMatrix(i,j) = -dualMatrix(i,j);


        if (!NTL::IsZero(dualMatrix(i,j) % matrix(j,j))) {
          BasInt gcd = NTL::GCD(dualMatrix(i,j), matrix(j,j));
          modulo *= matrix(j,j) / gcd;
          dualMatrix *= matrix(j,j) / gcd;
        }

        DivideRound (dualMatrix(i,j), matrix(j,j), dualMatrix(i,j));
      }

    }
  }

} // end namespace LatticeTester
