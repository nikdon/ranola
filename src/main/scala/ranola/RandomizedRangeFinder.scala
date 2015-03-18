package ranola


import breeze.linalg.{DenseMatrix, qr}
import breeze.stats.distributions.Rand


/**
 * Implementation of algorithms for computing an orthonormal matrix whose range approximates the range of M
 */
object RandomizedRangeFinder {

  /**
   * Algorithm 4.3 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   * 
   * @param M The input data matrix
   * @param sketchSize Size of the matrix to return
   * @param nIter Number of power iterations used to stabilize the result
   * @return A size-by-size projection matrix Q
   */
  def powerIteration(M: DenseMatrix[Double], sketchSize: Int, nIter: Int): DenseMatrix[Double] = {
    val R = DenseMatrix.rand(M.cols, sketchSize, rand = Rand.gaussian)
    val Y = M * R

    var a = 0
    while (a < nIter) {
      Y := M * (M.t * Y)
      a += 1
    }

    val q = qr.reduced.justQ(Y)
    q
  }

  /** Default parameter @param nIter = 10 */
  def powerIteration(M: DenseMatrix[Double], sketchSize: Int): DenseMatrix[Double] = {
    powerIteration(M, sketchSize, nIter = 10)
  }

  /** Default parameters: @param sketchSize = 10, @param nIter = 10 */
  def powerIteration(M: DenseMatrix[Double]): DenseMatrix[Double] = {
    powerIteration(M, sketchSize = 10, nIter = 10)
  }

  /**
   * Algorithm 4.1 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   * 
   * @param M The input data matrix
   * @param sketchSize Size of the matrix to return
   * @return A size-by-size projection matrix Q
   */
  def generic(M: DenseMatrix[Double], sketchSize: Int): DenseMatrix[Double] = {
    powerIteration(M, sketchSize, 0)
  }

  /** Default parameter @param sketchSize = 10 */
  def generic(M: DenseMatrix[Double]): DenseMatrix[Double] = {
    generic(M, 10)
  }
  
}
