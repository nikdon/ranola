package ranola


import breeze.linalg._
import breeze.math.Complex
import breeze.numerics.constants.Pi
import breeze.numerics.{exp, sqrt}
import breeze.signal.fourierTr
import breeze.stats.distributions.Rand


/**
 * Implementation of algorithms for computing an orthonormal matrix whose range approximates the range of M
 */
object RandomizedRangeFinder {

  /////////////////////////////////////////////////////////
  // Implementation of algorithms
  /////////////////////////////////////////////////////////

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
   * Algorithm 4.2 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param M The input data matrix
   * @param nRandVec Number of random vectors
   * @param tolerance Tolerance
   * @param maxIter Max number of iterations
   * @return
   */
  def adaptive(M: DenseMatrix[Double], nRandVec: Int, tolerance: Double, maxIter: Int): DenseMatrix[Double] = {
    require(nRandVec > 1, "Number of random vectors should be greater than 1")

    val (m, n) = (M.rows, M.cols)
    val threshold = tolerance / (10 * sqrt(2 / Math.PI))

    var W = new Array[DenseVector[Double]](nRandVec)
    var Y = new Array[DenseVector[Double]](nRandVec)
    var Q = DenseMatrix.create[Double](m, 0, new Array(0))

    // Helpers
    var i, j, k = 0
    val y = DenseVector.zeros[Double](m)
    val q = DenseVector.zeros[Double](m)
    val w = DenseVector.zeros[Double](n)
    val tmp = DenseVector.zeros[Double](m)

    // Fill W and Y arrays
    i = 0
    while (i < nRandVec) {
      W(i) = DenseVector.rand(n, rand = Rand.gaussian)
      Y(i) = M * W(i)
      i += 1
    }

    // Main loop
    i = 0
    j = 0
    while (maxNorm(Y) > threshold && i < maxIter) {
      j += 1

      y := Y(j) - Q * (Q.t * Y(j))
      q := (y / norm(y))

      if (i == 0) Q = q.toDenseMatrix.t
      else Q = DenseMatrix.horzcat(Q, q.asDenseMatrix.t)

      w := DenseVector.rand(n, rand = Rand.gaussian)
      tmp := M * w
      y := tmp - Q * (Q.t * tmp)

      W :+= w
      Y :+= y

      // Overwrite Y
      k = j + 1
      while (k < j + nRandVec) {
        Y(i) = Y(i) - (q.t * Y(i)) * q
        k += 1
      }

      i += 1
    }

    Q
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
   * Algorithm 4.5 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param M The input data matrix
   * @param sketchSize Size of the matrix to return
   * @return A size-by-size projection matrix Q
   */
  def fastGeneric(M: DenseMatrix[Double], sketchSize: Int): DenseMatrix[Double] = {
    val (m, n) = (M.rows, M.cols)
    val D = ??? // Diagonal matrix with complex var uniformly distributed on complex unit circle
    val F = ??? // DFT matrix
    val R = ??? // Matrix with random columns of the I
    val SRFT = ??? // sqrt(n / sketchSize) * D * F * R
    val Y = ??? // M * SRFT   // See $3.3 in [[http://www.cs.yale.edu/homes/el327/papers/approximationOfMatrices.pdf]]
    val q = ??? // qr.reduced.justQ(Y)
    ??? // q
  }

  def D(n: Int): DenseMatrix[Complex] = {
    diag(DenseVector(Array.tabulate(n)(a => exp(Complex.i * 2 * Pi * util.Random.nextGaussian()))))
  }

  def F(n: Int): DenseMatrix[Complex] = {
    val x = DenseMatrix.zeros[Double](n, n)

    for {
      i <- 0 until n
      j <- 0 until n
    } x(i, j) = ???

    fourierTr(x)
  }

  def R(n: Int, sketchSize: Int): DenseMatrix[Double] = {
    require(sketchSize <= n)

    val data = Array.tabulate(n)(a => a)
    val samples = shuffle(data, sketchSize).toSeq
    val I = DenseMatrix.eye[Double](n)

    I(::, samples).toDenseMatrix
  }

  // Fisher-Yates shuffle
  // See: [[http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle]]
  def shuffle[@specialized(Double) T](array: Array[T],
                                      sketchSize: Int,
                                      swap: (T, T) => (T, T) = { (a: T, b: T) => (b, a) }): Array[T] = {
    val random = new scala.util.Random

    for (n <- array.length - 1 to 0 by -1) {
      val k = random.nextInt(n + 1)
      val (a, b) = swap(array(k), array(n))
      array(k) = a
      array(n) = b
    }

    array.slice(0, sketchSize)
  }


  ///////////////////////////////////////
  // Helpers
  ///////////////////////////////////////

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

  /**
   * Maximum of DenseVector norms in the array
   *
   * @param y Array contained DenseVectors
   * @return Maximum norm of the given array
   */
  private def maxNorm(y: Array[DenseVector[Double]]): Double = {
    val len = y.length
    val norms = new Array[Double](len)
    var i = 0
    while (i < len) {
      norms(i) = norm(y(i))
      i += 1
    }
    max(norms)
  }

}
