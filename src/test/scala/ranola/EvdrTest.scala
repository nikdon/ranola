package ranola


import breeze.linalg.eigSym.EigSym
import breeze.linalg.{CSCMatrix, DenseMatrix, DenseVector, argsort}
import breeze.util.DoubleImplicits
import org.junit.runner.RunWith
import org.scalatest._
import org.scalatest.junit._
import ranola.TestHelpers._


@RunWith(classOf[JUnitRunner])
class EvdrTest extends FunSuite with Matchers with DoubleImplicits {

  val A_sparse = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
  val A_dense = A_sparse.toDense
  
  def checkResults(evals: DenseVector[Double], evect: DenseMatrix[Double]): Unit = {
    val idx = argsort(evals)

    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    idx.zipWithIndex.foreach { case (i, n) =>
      evals(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evect(::, i), eigVect(::, n), 1E-6)
    }
  }

  ///////////////////////
  // DENSE
  ///////////////////////

  test("EVDR with Generic Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.generic(A_dense, k = 3, overSamples = 2)
    checkResults(lambda, evs)
  }

  test("EVDR with Fast Generic Range Finder") {
    val EigSym(lambda, evs) = evdr.fastGeneric(A_dense, k = 3, overSamples = 0)
    checkResults(lambda, evs)
  }

  test("EVDR with Power Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.powerIteration(A_dense, k = 3, nIter = 5, overSamples = 2)
    checkResults(lambda, evs)
  }

  test("EVDR with Subspace Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.subspaceIteration(A_dense, k = 3, nIter = 5, overSamples = 2)
    checkResults(lambda, evs)
  }

  test("EVDR with Adaptive Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.adaptive(A_dense, k = 3, overSamples = 0, tol = 1E-3, maxIter = 3)
    checkResults(lambda, evs)
  }

  
  ///////////////////////
  // SPARSE
  ///////////////////////

  test("Sparse EVDR with Generic Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.generic(A_sparse, k = 3, overSamples = 2)
    checkResults(lambda, evs)
  }

  test("Sparse EVDR with Fast Generic Range Finder") {
    val EigSym(lambda, evs) = evdr.fastGeneric(A_sparse, k = 3, overSamples = 0)
    checkResults(lambda, evs)
  }

  test("Sparse EVDR with Power Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.powerIteration(A_sparse, k = 3, nIter = 5, overSamples = 2)
    checkResults(lambda, evs)
  }

  test("Sparse EVDR with Subspace Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.subspaceIteration(A_sparse, k = 3, nIter = 5, overSamples = 2)
    checkResults(lambda, evs)
  }

  test("Sparse EVDR with Adaptive Iteration Randomized Range Finder") {
    val EigSym(lambda, evs) = evdr.adaptive(A_sparse, k = 3, overSamples = 0, tol = 1E-3, maxIter = 3)
    checkResults(lambda, evs)
  }
}
