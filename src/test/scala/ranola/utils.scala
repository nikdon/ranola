package ranola


import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.util.DoubleImplicits
import org.scalatest._

object utils extends FunSuite with Matchers with DoubleImplicits {

  def vectorsNearlyEqual(A: DenseVector[Double], B: DenseVector[Double], threshold: Double = 1E-6) {
    for(i <- 0 until A.length)
      A(i) should be (B(i) +- threshold)
  }

  def matricesNearlyEqual(A: DenseMatrix[Double], B: DenseMatrix[Double], threshold: Double = 1E-6) {
    for(i <- 0 until A.rows; j <- 0 until A.cols)
      A(i,j) should be (B(i, j) +- threshold)
  }

}
