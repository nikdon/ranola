package ranola


import breeze.linalg._
import breeze.linalg.svd._
import breeze.numerics.abs
import breeze.util.DoubleImplicits
import org.junit.runner.RunWith
import org.scalatest._
import org.scalatest.junit._
import ranola.TestHelpers._


@RunWith(classOf[JUnitRunner])
class SvdrTest extends FunSuite with Matchers with DoubleImplicits {

  ///////////////////////
  // DENSE
  ///////////////////////

  test("svd and svdr.generic singular values are equal") {
    val a = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val SVD(u, s, v) = svd.reduced(m)
      val SVD(ur, sr, vr) = svdr.generic(m, k = m.rows min m.cols, overSamples = 1)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
      matricesNearlyEqual(m, ur * diag(sr) * vr)
    }
  }

  test("svd and svdr.fastGeneric singular values are equal") {
    val a = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val SVD(u, s, v) = svd.reduced(m)
      val SVD(ur, sr, vr) = svdr.fastGeneric(m, k = m.rows min m.cols, overSamples = 0)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
      matricesNearlyEqual(m, ur * diag(sr) * vr)
    }
  }

  test("svd and svdr.powerIteration singular values are equal") {
    val a = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val SVD(u, s, v) = svd.reduced(m)
      val SVD(ur, sr, vr) = svdr.powerIteration(m, k = m.rows min m.cols, nIter = 5, overSamples = 1)


      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
      matricesNearlyEqual(m, ur * diag(sr) * vr)
    }
  }

  test("svd and svdr.adaptive singular values are equal A[m x n], m > n") {
    // Note: Size of Q depends on maxIter and overSamples,
    // so these should be true:
    // (maxIter + overSamples) <= A.cols   ==> for convergence
    // (maxIter + overSamples) > k         ==> for randomized decomposition

    val a = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    val SVD(u, s, v) = svd.reduced(a)
    val SVD(ur, sr, vr) = svdr.adaptive(a, k = a.rows min a.cols, tol = 1E-3, maxIter = 5, overSamples = 0)


    vectorsNearlyEqual(s, sr)
    matricesNearlyEqual(abs(u), abs(ur))
    matricesNearlyEqual(abs(v), abs(vr))
    matricesNearlyEqual(a, ur * diag(sr) * vr)

  }

  test("svd and svdr.adaptive singular values are equal for A[m x n], m < n") {
    // Note: Size of Q depends on maxIter and overSamples,
    // so these should be true:
    // (maxIter + overSamples) <= A.cols   ==> for convergence
    // (maxIter + overSamples) > k         ==> for randomized decomposition

    val a = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    ).t

    val SVD(u, s, v) = svd.reduced(a)
    val SVD(ur, sr, vr) = svdr.adaptive(a, k = a.rows min a.cols, tol = 1E-3, maxIter = 3, overSamples = 0)


    vectorsNearlyEqual(s, sr)
    matricesNearlyEqual(abs(u), abs(ur))
    matricesNearlyEqual(abs(v), abs(vr))
    matricesNearlyEqual(a, ur * diag(sr) * vr)

  }

  test("svd and svdr.subspaceIteration singular values are equal") {
    val a = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val SVD(u, s, v) = svd.reduced(m)
      val SVD(ur, sr, vr) = svdr.subspaceIteration(m, k = m.rows min m.cols, nIter = 5, overSamples = 5)


      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
      matricesNearlyEqual(m, ur * diag(sr) * vr)
    }
  }


  ///////////////////////
  // SPARSE
  ///////////////////////


  test("Sparse svd and svdr.generic singular values are equal") {
    val a = CSCMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val k = 2
      val SVD(u, s, v) = svd(m, k, 1E-6)
      val SVD(ur, sr, vr) = svdr.generic(m, k, overSamples = 1)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
    }
  }

  test("Sparse svd and svdr.fastGeneric singular values are equal") {
    val a = CSCMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val k = 2
      val SVD(u, s, v) = svd(m, k, 1E-6)
      val SVD(ur, sr, vr) = svdr.fastGeneric(m, k, overSamples = 1)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
    }
  }

  test("Sparse svd and svdr.powerIteration singular values are equal") {
    val a = CSCMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val k = 2
      val SVD(u, s, v) = svd(m, k, 1E-6)
      val SVD(ur, sr, vr) = svdr.powerIteration(m, k, nIter = 5, overSamples = 1)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
    }
  }

  test("Sparse svd and svdr.adaptive singular values are equal A[m x n], m > n") {
    val a = CSCMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    val k = 2
    val SVD(u, s, v) = svd(a, k, 1E-6)
    val SVD(ur, sr, vr) = svdr.adaptive(a, k, tol = 1E-3, maxIter = 3, overSamples = 0)

    vectorsNearlyEqual(s, sr)
    matricesNearlyEqual(abs(u), abs(ur))
    matricesNearlyEqual(abs(v), abs(vr))
  }

  test("Sparse svd and svdr.adaptive singular values are equal for A[m x n], m < n") {
    val a = CSCMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    ).t

    val k = 2
    val SVD(u, s, v) = svd(a, k, 1E-6)
    val SVD(ur, sr, vr) = svdr.adaptive(a, k, tol = 1E-3, maxIter = 3, overSamples = 0)

    vectorsNearlyEqual(s, sr)
    matricesNearlyEqual(abs(u), abs(ur))
    matricesNearlyEqual(abs(v), abs(vr))
  }

  test("Sparse svd and svdr.subspaceIteration singular values are equal") {
    val a = CSCMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    for (m <- List(a, a.t)) {
      val k = 2
      val SVD(u, s, v) = svd(m, k, 1E-6)
      val SVD(ur, sr, vr) = svdr.subspaceIteration(m, k, nIter = 5, overSamples = 5)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
    }
  }
}
