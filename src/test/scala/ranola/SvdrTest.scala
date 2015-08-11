package ranola


import breeze.linalg._
import breeze.linalg.svd.SVD
import breeze.numerics.abs
import breeze.util.DoubleImplicits
import org.scalatest._


class SvdrTest extends FunSuite with Matchers with DoubleImplicits {

  val A_dense = DenseMatrix(
    (2.0, 4.0, 0.0),
    (1.0, 3.0, 4.0),
    (5.0, 0.0, 0.9),
    (3.0, 5.0, 0.5),
    (7.5, 1.0, 6.0),
    (0.0, 7.0, 0.0)
  )

  def vectorsNearlyEqual(A: DenseVector[Double], B: DenseVector[Double], threshold: Double = 1E-6): Unit = {
    for (i <- 0 until A.length)
      A(i) should be(B(i) +- threshold)
  }

  def matricesNearlyEqual(A: DenseMatrix[Double], B: DenseMatrix[Double], threshold: Double = 1E-6): Unit = {
    for (i <- 0 until A.rows; j <- 0 until A.cols)
      A(i, j) should be(B(i, j) +- threshold)
  }

  def checkResults(evals: DenseVector[Double], evect: DenseMatrix[Double]): Unit = {
    val idx = argsort(evals)

    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    idx.zipWithIndex.foreach { case (i, n) =>
      evals(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evect(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("SVD and SVDR with Power Iteration Randomized Range Finder") {
    for (m <- List(A_dense, A_dense.t)) {
      val SVD(u, s, v) = svd.reduced(m)
      val SVD(ur, sr, vr) = SVDR.viaPowerIteration(m, k = m.rows min m.cols, nIter = 5, overSamples = 1)


      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
      matricesNearlyEqual(m, ur * diag(sr) * vr)
    }
  }
}