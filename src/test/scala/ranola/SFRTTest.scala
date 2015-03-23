package ranola


import breeze.util.DoubleImplicits
import org.junit.runner.RunWith
import org.scalatest._
import org.scalatest.junit._


import ranola.RandomizedRangeFinder._


@RunWith(classOf[JUnitRunner])
class SFRTTest extends FunSuite with Matchers with DoubleImplicits {

  test("F Matrix") {
    val matF = F(3)
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
}
