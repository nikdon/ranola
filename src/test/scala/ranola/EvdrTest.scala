package ranola


import breeze.linalg.eigSym.EigSym
import breeze.linalg.{CSCMatrix, DenseMatrix, DenseVector, argsort}
import breeze.util.DoubleImplicits
import org.scalatest._
import ranola.EVDR.Implicits._


class EvdrTest extends FunSuite with Matchers with DoubleImplicits {

  val A_sparse = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
  val A_dense = A_sparse.toDense

  def vectorsNearlyEqual(A: DenseVector[Double], B: DenseVector[Double], threshold: Double = 1E-6): Unit = {
    for (i <- 0 until A.length)
      A(i) should be(B(i) +- threshold)
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

  test("EVDR via Generic Randomized Range Finder") {
    val EigSym(lambda, evs) = EVDR.generic(A_dense, k = 3, nOverSamples = 2)
    checkResults(lambda, evs)
  }

  test("EVDR via Power Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = EVDR.viaPowerIteration(A_dense, k = 3, nIter = 5, nOverSamples = 2)
    checkResults(lambda, evs)
  }
}