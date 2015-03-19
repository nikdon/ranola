package ranola


import breeze.linalg.{DenseMatrix, qr}
import breeze.stats.distributions.Rand


/**
 * Implementation of algorithms for computing an orthonormal matrix whose range approximates the range of M
 */
object RandomizedRangeFinder {

  private val DEFAULT_N_ITER = 10

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
    val R = drawRandomMatrix(M, sketchSize)
    val Y = M * R

    var i = 0
    while (i < nIter) {
      Y := M * (M.t * Y)
      i += 1
    }

    val q = qr.reduced.justQ(Y)
    q
  }

  def powerIteration(M: DenseMatrix[Double], sketchSize: Int): DenseMatrix[Double] = {
    powerIteration(M, sketchSize, nIter = DEFAULT_N_ITER)
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

  /**
   * Algorithm 4.4 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param M The input data matrix
   * @param sketchSize Size of the matrix to return
   * @return A size-by-size projection matrix Q
   */
  def subspaceIteration(M: DenseMatrix[Double], sketchSize: Int, nIter: Int): DenseMatrix[Double] = {
    val R = drawRandomMatrix(M, sketchSize)
    val Y = M * R
    val q = qr.reduced.justQ(Y)

    var i = 0
    while (i < nIter) {
      Y := M.t * q
      q := qr.reduced.justQ(Y)
      Y := M * q
      q := qr.reduced.justQ(Y)
      i += 1
    }

    q
  }

  def subspaceIteration(M: DenseMatrix[Double], sketchSize: Int): DenseMatrix[Double] = {
    subspaceIteration(M, sketchSize, nIter = DEFAULT_N_ITER)
  }

  /**
   * Draw random matrix
   *
   * @param M Matrix
   * @param sketchSize proposed sketch size that is second dimension of random matrix
   * @return Random matrix
   */
  private def drawRandomMatrix(M: DenseMatrix[Double], sketchSize: Int): DenseMatrix[Double] = {
    val n = checkAndGetSketchSize(M, sketchSize)
    val l = M.cols
    DenseMatrix.rand(l, n, rand = Rand.gaussian)
  }

  /**
   * Checks the proposed size of the sketch.
   * @param M Matrix to be sampled
   * @param s proposed sketch size
   * @return Correct sketch size
   */
  private def checkAndGetSketchSize(M: DenseMatrix[Double], s: Int): Int = {
    if (s > M.rows) M.rows
    else s
  }

}
