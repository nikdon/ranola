package ranola

import breeze.generic.UFunc
import breeze.linalg._
import breeze.linalg.eig.Eig
import breeze.linalg.eigSym.{DenseEigSym, EigSym}
import breeze.numerics._
import ranola.utils._


/**
 * Approximate truncated randomized EVD
 */
object evdr extends UFunc {

  implicit object EVDR_DM_Impl2 extends Impl2[DenseMatrix[Double], Int, DenseEigSym] {
    def apply(M: DenseMatrix[Double], s: Int): DenseEigSym =
      doEigSymDouble(M, s, nOversamples = 10, nIter = 0)
  }

  implicit object EVDR_DM_Impl3 extends Impl3[DenseMatrix[Double], Int, Int, DenseEigSym] {
    def apply(M: DenseMatrix[Double], s: Int, nOversamples: Int): DenseEigSym =
      doEigSymDouble(M, s, nOversamples, nIter = 0)
  }

  implicit object EVDR_DM_Impl4 extends Impl4[DenseMatrix[Double], Int, Int, Int, DenseEigSym] {
    def apply(M: DenseMatrix[Double], s: Int, nOversamples: Int, nIter: Int): DenseEigSym =
      doEigSymDouble(M, s, nOversamples, nIter)
  }

  /**
   * Computes an approximate truncated randomized EVD. Fast on large matrices.
   *
   * @param M Matrix to decompose
   * @param s Number of columns in orthonormal matrix (sketch size)
   * @param nOversamples Additional number of random vectors to sample the range of M so as
   *                     to ensure proper conditioning. The total number of random vectors
   *                     used to find the range of M is [s + nOversamples]
   * @param nIter Number of power iterations (can be used to deal with very noisy problems)
   * @return The eigenvalue decomposition (EVD) with the eigenvalues and the eigenvectors
   *
   * ==References==
   *
   * Finding structure with randomness: Stochastic algorithms for constructing
   * approximate matrix decompositions
   * Halko, et al., 2009 [[http://arxiv.org/abs/arXiv:0909.4061]]
   */
  private def doEigSymDouble(M: DenseMatrix[Double],
                             s: Int,
                             nOversamples: Int = 10,
                             nIter: Int = 0): DenseEigSym = {

    require(s <= (M.rows min M.cols), "Number of columns in orthonormal matrix should be less than min(M.rows, M.cols)")
    require(s >= 1, "Sketch size should be greater than 1")

    val nRandom = s + nOversamples

    val Q = RandomizedRangeFinder.powerIteration(M, nRandom, nIter)

    val b = Q.t * (M * Q)

    val Eig(w, _, v) = eig(b)

    val _u = Q * v

    val u = utils.flipSigns(_u)

    EigSym(w, u)
  }


}
