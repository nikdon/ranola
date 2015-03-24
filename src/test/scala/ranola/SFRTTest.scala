package ranola


import breeze.linalg.{qr, DenseMatrix}
import breeze.math.Complex
import breeze.numerics._
import breeze.util.DoubleImplicits
import org.junit.runner.RunWith
import org.scalatest._
import org.scalatest.junit._


import ranola.RandomizedRangeFinder._


@RunWith(classOf[JUnitRunner])
class SFRTTest extends FunSuite with Matchers with DoubleImplicits {

  test("F Matrix") {
    val matF = F(4)
    println(matF)
  }

  test("D Matrix") {
    val matD = D(3)
    println(matD)
  }

  test("R Matrix") {
    val matR = R(5, 3)
    println(matR)
  }

  test("SRFT Matrix") {
    val n = 3
    val s = 3

    val d = D(n)
    val f = F(n)
    val r = R (n , s)
    val srft = d * f * r * Complex(sqrt(n / s), 0)

    println(srft)
  }
}
