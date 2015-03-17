package ranola


import breeze.linalg.{qr, DenseMatrix}
import breeze.stats.distributions.Rand


/**
 * Implementation of algorithms for computing an orthonormal matrix whose range approximates the range of M
 */
object RandomizedRangeFinder {

  /**
   * Algorithm 4.1 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param M The input data matrix
   * @param size Size of the matrix to return
   * @return A size-by-size projection matrix Q
   */
  def generic(M: DenseMatrix[Double], size: Int) = {
    powerIteration(M, size, 0)
  }

  /**
   * Algorithm 4.3 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param M The input data matrix
   * @param size Size of the matrix to return
   * @param nIter Number of power iterations used to stabilize the result
   * @return A size-by-size projection matrix Q
   */
  def powerIteration(M: DenseMatrix[Double],
                     size: Int,
                     nIter: Int): DenseMatrix[Double] = {
    val R = DenseMatrix.rand(M.cols, size, rand = Rand.gaussian)
    val Y = M * R
    for (a <- 0 until nIter) Y := M * (M.t * Y)
    val q = qr.reduced.justQ(Y)
    q
  }

}
