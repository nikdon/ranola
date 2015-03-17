package ranola


import breeze.linalg.svd.SVD
import breeze.numerics.abs
import org.scalatest._
import org.scalatest.junit._
import breeze.linalg._
import breeze.util.DoubleImplicits
import org.junit.runner.RunWith

import ranola.utils._


@RunWith(classOf[JUnitRunner])
class SvdrTest extends FunSuite with Matchers with DoubleImplicits {

  test("svd and svdr singular values are equal") {
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
      val SVD(ur, sr, vr) = svdr(m, m.rows min m.cols)

      vectorsNearlyEqual(s, sr)
      matricesNearlyEqual(abs(u), abs(ur))
      matricesNearlyEqual(abs(v), abs(vr))
    }
  }

  test("svdr A[m, n], m < n") {
    val m = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    ).t

    val SVD(u, sr, vt) = svdr(m, m.rows min m.cols)

    val reM = u * diag(sr) * vt
    matricesNearlyEqual(reM, m)
  }

  test("svdr A[m, n], m > n") {
    val m = DenseMatrix(
      (2.0, 4.0, 0.0),
      (1.0, 3.0, 4.0),
      (5.0, 0.0, 0.9),
      (3.0, 5.0, 0.5),
      (7.5, 1.0, 6.0),
      (0.0, 7.0, 0.0)
    )

    val SVD(u, sr, vt) = svdr(m, m.rows min m.cols)

    val reM = u * diag(sr) * vt
    matricesNearlyEqual(reM, m)
  }


}
