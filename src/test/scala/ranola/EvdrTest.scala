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

  ///////////////////////
  // DENSE
  ///////////////////////

  test("EVDR with Generic Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.generic(A, k = 3, overSamples = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("EVDR with Fast Generic Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.fastGeneric(A, k = 3, overSamples = 0)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("EVDR with Power Iteration Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.powerIteration(A, k = 3, nIter = 5, overSamples = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("EVDR with Subspace Iteration Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.subspaceIteration(A, k = 3, nIter = 5, overSamples = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("EVDR with Adaptive Iteration Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.adaptive(A, k = 3, tol = 1E-3, maxIter = 3, overSamples = 0)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  ///////////////////////
  // SPARSE
  ///////////////////////

  test("Sparse EVDR with Generic Randomized Range Finder") {
    val A = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.generic(A, k = 3, overSamples = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("Sparse EVDR with Fast Generic Range Finder") {
    val A = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.fastGeneric(A, k = 3, overSamples = 0)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("Sparse EVDR with Power Iteration Randomized Range Finder") {
    val A = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.powerIteration(A, k = 3, nIter = 5, overSamples = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("Sparse EVDR with Subspace Iteration Randomized Range Finder") {
    val A = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.subspaceIteration(A, k = 3, nIter = 5, overSamples = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }

  test("Sparse EVDR with Adaptive Iteration Randomized Range Finder") {
    val A = CSCMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0, 25.0, 82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr.adaptive(A, k = 3, tol = 1E-3, maxIter = 3, overSamples = 0)

    val idx = argsort(lambda)

    idx.zipWithIndex.foreach { case (i, n) =>
      lambda(i) should be(eigVals(n) +- 1E-6)
      vectorsNearlyEqual(evs(::, i), eigVect(::, n), 1E-6)
    }
  }
}
