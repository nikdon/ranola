package ranola

import breeze.linalg.qr.QR
import breeze.linalg.{qr, DenseVector, DenseMatrix}
import breeze.stats.distributions.Rand


/**
 * Trait contains common matrix operations
 *
 * @tparam N The type og values inside matrix and vector containers
 * @tparam M The matrix container
 * @tparam V The vector container
 */
trait MatrixOps[N, M[_], V[_]] {

  def absM(A: M[N]): M[N]
  def getCols(A: M[N]): Int
  def getRows(A: M[N]): Int

  /** Matrix-matrix multiplication */
  def mltMV(A: M[N], B: V[N]): M[N] // A * B

  /** Matrix-vector multiplication */
  def mltMM(A: M[N], B: M[N]): M[N] // A * B

  /** Matrix transpose */
  def t(A: M[N]): M[N]

  /** Create a new matrix filled with random values */
  def drawRandomMatrix(A: M[N], sketchSize: Int): M[N]

  /** Create a new matrix filled with zeros */
  def drawZerosMatrix(nRow: Int, nCol: Int): M[N]

  /** Check size of the given matrix and return an appropriate sketch size */
  protected [this] def checkAndGetSketchSize(A: M[N], s: Int): Int = if (s > getRows(A)) getRows(A) else s

  /** QR decomposition of the given matrix A */
  def QR(A: M[N]): (M[N], M[N])

}


object MatrixOps {

  //////////////////////////////////
  // Add new implementations here //
  //////////////////////////////////

  implicit val BreezeDoubleMatrixOps = new MatrixOps[Double, breeze.linalg.DenseMatrix, breeze.linalg.DenseVector] {
    override def absM(A: DenseMatrix[Double]): DenseMatrix[Double] = breeze.numerics.abs(A)
    override def getCols(A: DenseMatrix[Double]) = A.cols
    override def getRows(A: DenseMatrix[Double]) = A.rows

    override def mltMV(A: DenseMatrix[Double], B: DenseVector[Double]): DenseMatrix[Double] = A * B.toDenseMatrix
    override def mltMM(A: DenseMatrix[Double], B: DenseMatrix[Double]): DenseMatrix[Double] = A * B
    override def t(A: DenseMatrix[Double]): DenseMatrix[Double] = A.t

    override def drawRandomMatrix(A: DenseMatrix[Double], sketchSize: Int) = {
      val n = checkAndGetSketchSize(A, sketchSize)
      val l = getCols(A)
      DenseMatrix.rand(l, n, rand = Rand.gaussian)
    }

    override def drawZerosMatrix(nRow: Int, nCol: Int) = DenseMatrix.zeros[Double](nRow, nCol)

    override def QR(A: DenseMatrix[Double]) = {
      val QR(_q, _r) = qr.reduced(A)
      (_q, _r)
    }
  }

}
