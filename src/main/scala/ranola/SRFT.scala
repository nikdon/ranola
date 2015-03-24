package ranola


import breeze.linalg.{*, DenseMatrix, DenseVector, diag}
import breeze.math.Complex
import breeze.numerics._
import breeze.numerics.constants.Pi
import breeze.signal.fourierTr


/**
 * A Subsampled Random Fourier Transform. A random matrix with structure.
 */
object SRFT {

  def apply(n: Int, sketchSize: Int): DenseMatrix[Complex] = {
    val d = SRFT.D(n)
    val f = SRFT.F(n)
    val r = SRFT.R(n, sketchSize)
    val srft = (d * (f * r)) * Complex(sqrt(n / sketchSize), 0)
    srft
  }

  /**
   * Constructor of a diagonal matrix with complex var uniformly distributed on complex unit circle
   *
   * @param n Size
   * @return A n-by-n diagonal matrix
   */
  def D(n: Int): DenseMatrix[Complex] = {
    diag(DenseVector(Array.tabulate(n)(a => exp(Complex.i * 2 * Pi * util.Random.nextGaussian()))))
  }

  /**
   * Constructor of a DFT matrix
   *
   * @param n Size
   * @return A n-by-n matrix
   */
  def F(n: Int): DenseMatrix[Complex] = {
    val I = DenseMatrix.eye[Complex](n)
    I(*, ::).map(dv => fourierTr(dv))
  }

  /**
   * Constructor of a matrix with random columns of the I
   * @param n Size
   * @param sketchSize Sketch size
   * @return A n-by-sketchSize matrix
   */
  def R(n: Int, sketchSize: Int): DenseMatrix[Complex] = {
    require(sketchSize <= n)

    val I = DenseMatrix.eye[Complex](n)
    val data = Array.tabulate(n)(a => a)
    val samples = shuffle(data, sketchSize).toSeq

    I(::, samples).toDenseMatrix
  }

  /**
   * Fisher-Yates shuffle. See [[http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle]]
   *
   * @param array Array to shuffle
   * @param sketchSize Size of the sketch
   * @param swap Swap function
   * @tparam T Any
   * @return Shuffled array sliced from 0 up to sketchSize
   */
  private def shuffle[T](array: Array[T],
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

}
