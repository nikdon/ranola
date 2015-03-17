package ranola

import breeze.generic.UFunc
import breeze.linalg._
import breeze.linalg.svd.{DenseSVD, SVD}
import breeze.numerics.{abs, signum}
import breeze.stats.distributions.Rand


/**
 * Approximate truncated randomized SVD
 */
object svdr extends UFunc {

  implicit object SvdR_DM_Impl2 extends Impl2[DenseMatrix[Double], Int, DenseSVD] {
    def apply(M: DenseMatrix[Double], k: Int): DenseSVD =
      doSVDR_Double(M, k, nOversamples = 10, nIter = 0)
  }

  implicit object SvdR_DM_Impl3 extends Impl3[DenseMatrix[Double], Int, Int, DenseSVD] {
    def apply(M: DenseMatrix[Double], k: Int, nOversamples: Int): DenseSVD =
      doSVDR_Double(M, k, nOversamples, nIter = 0)
  }

  implicit object SvdR_DM_Impl4 extends Impl4[DenseMatrix[Double], Int, Int, Int, DenseSVD] {
    def apply(M: DenseMatrix[Double], k: Int, nOversamples: Int, nIter: Int): DenseSVD =
      doSVDR_Double(M, k, nOversamples, nIter)
  }

  /**
   * Computes an approximate truncated randomized SVD. Fast on large matrices with limited
   * number of singular values and vectors that are necessary to extract
   *
   * @param M Matrix to decompose
   * @param k Number of singular values and vectors to extract
   * @param nOversamples Additional number of random vectors to sample the range of M so as
   *                     to ensure proper conditioning. The total number of random vectors
   *                     used to find the range of M is [k + nOversamples]
   * @param nIter Number of power iterations (can be used to deal with very noisy problems)
   * @return The singular value decomposition (SVD) with the singular values,
   *         the left and right singular vectors
   *
   * ==References==
   *
   * Finding structure with randomness: Stochastic algorithms for constructing
   * approximate matrix decompositions
   * Halko, et al., 2009 [[http://arxiv.org/abs/arXiv:0909.4061]]
   */
  private def doSVDR_Double(M: DenseMatrix[Double],
                            k: Int,
                            nOversamples: Int = 10,
                            nIter: Int = 0): DenseSVD = {

    require(k <= (M.rows min M.cols), "Number of singular values should be less than min(M.rows, M.cols)")

    val nRandom = k + nOversamples

    val Q = RandomizedRangeFinder.powerIteration(M, nRandom, nIter)

    val b = Q.t * M

    val SVD(w2, _s, _v) = svd.reduced(b)

    val _u = Q * w2

    val (u, v) = utils.flipSVDSigns(_u, _v)

    SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
  }
}
