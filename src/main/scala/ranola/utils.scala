package ranola


import breeze.linalg.{DenseMatrix, argmax}
import breeze.numerics.{abs, signum}


object utils {

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive.
   * Applicable for EVD.
   *
   * @param u eigenvectors
   * @return eigenvectors with resolved sign ambiguity
   */
  def flipEVDSigns(u: DenseMatrix[Double]): DenseMatrix[Double] = {
    val abs_u = abs(u)
    val len = u.cols
    var i = 0
    while (i < len) {
      u(::, i) :*= signum(u(argmax(abs_u(::, i)), i))
      i += 1
    }

    u
  }

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive.
   * Applicable for SVD.
   *
   * @param u left singular vectors
   * @param v right singular vectors
   * @return left and right singular vectors with resolved sign ambiguity
   */
  def flipSVDSigns(u: DenseMatrix[Double],
                   v: DenseMatrix[Double]): (DenseMatrix[Double], DenseMatrix[Double]) = {
    val abs_u = abs(u)
    val len = u.cols
    var i = 0
    while (i < len) {
      val sign = signum(u(argmax(abs_u(::, i)), i))
      u(::, i) :*= sign
      v(i, ::) :*= sign
      i += 1
    }

    (u, v)
  }

}
