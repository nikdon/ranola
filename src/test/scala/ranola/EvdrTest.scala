package ranola

import org.scalatest._
import org.scalatest.junit._
import breeze.linalg.eigSym.EigSym
import breeze.linalg.{argsort, DenseVector, DenseMatrix}
import breeze.util.DoubleImplicits
import org.junit.runner.RunWith

import ranola.TestHelpers._


@RunWith(classOf[JUnitRunner])
class EvdrTest extends FunSuite with Matchers with DoubleImplicits {

  test("EVDR with Generic Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0,25.0,82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr(A, RandomizedRangeFinder.generic(A, sketchSize = 2), k = 2)

    val idx = argsort(lambda)

    idx.zipWithIndex.map{ i =>
      lambda(i._1) should be (eigVals(i._2) +- 1E-6)
      vectorsNearlyEqual(evs(::, i._1), eigVect(::, i._2), 1E-6)
    }
  }

  test("EVDR with Power Iteration Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0,25.0,82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr(A, RandomizedRangeFinder.powerIteration(A, sketchSize = 3), k = 3)

    val idx = argsort(lambda)

    idx.zipWithIndex.map{ i =>
      lambda(i._1) should be (eigVals(i._2) +- 1E-6)
      vectorsNearlyEqual(evs(::, i._1), eigVect(::, i._2), 1E-6)
    }
  }

  test("EVDR with Subspace Iteration Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0,25.0,82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr(A, RandomizedRangeFinder.subspaceIteration(A, sketchSize = 3), k = 3)

    val idx = argsort(lambda)

    idx.zipWithIndex.map{ i =>
      lambda(i._1) should be (eigVals(i._2) +- 1E-6)
      vectorsNearlyEqual(evs(::, i._1), eigVect(::, i._2), 1E-6)
    }
  }

  test("EVDR with Adaptive Iteration Randomized Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0,25.0,82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr(A, RandomizedRangeFinder.adaptive(A, nRandVec = 2, tolerance = 0.01, maxIter = 3), k = 3)

    val idx = argsort(lambda)

    idx.zipWithIndex.map{ i =>
      lambda(i._1) should be (eigVals(i._2) +- 1E-6)
      vectorsNearlyEqual(evs(::, i._1), eigVect(::, i._2), 1E-6)
    }
  }

  test("EVDR with Fast Generic Range Finder") {
    val A = DenseMatrix((9.0, 0.0, 0.0), (0.0, 82.0, 0.0), (0.0, 0.0, 25.0))
    val eigVals = DenseVector(9.0,25.0,82.0)
    val eigVect = DenseMatrix((1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0))

    val EigSym(lambda, evs) = evdr(A, RandomizedRangeFinder.fastGeneric(A, sketchSize = 3), k = 3)

    val idx = argsort(lambda)

    idx.zipWithIndex.map{ i =>
      lambda(i._1) should be (eigVals(i._2) +- 1E-6)
      vectorsNearlyEqual(evs(::, i._1), eigVect(::, i._2), 1E-6)
    }
  }
}
