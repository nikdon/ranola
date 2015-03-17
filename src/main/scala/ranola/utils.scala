package ranola


import breeze.linalg.{argmax, DenseMatrix}
import breeze.numerics._


object utils {

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive
   *
   * @param u eigenvectors
   * @return eigenvectors with resolved sign ambiguity
   */
  def flipSigns(u: DenseMatrix[Double]): DenseMatrix[Double] = {
    val abs_u = abs(u)
    val max_abs_cols = (0 until u.cols).map(c => argmax(abs_u(::, c)))
    val signs = max_abs_cols.zipWithIndex.map(e => signum(u(e._1, e._2)))
    signs.zipWithIndex.foreach(s => {
      u(::, s._2) :*= s._1
    })
    u
  }

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive
   *
   * @param u left singular vectors
   * @param v right singular vectors
   * @return left and right singular vectors with resolved sign ambiguity
   */
  def flipSVDSigns(u: DenseMatrix[Double],
                   v: DenseMatrix[Double]): (DenseMatrix[Double], DenseMatrix[Double]) = {
    val abs_u = abs(u)
    val max_abs_cols = (0 until u.cols).map(c => argmax(abs_u(::, c)))
    val signs = max_abs_cols.zipWithIndex.map(e => signum(u(e._1, e._2)))
    signs.zipWithIndex.foreach(s => {
      u(::, s._2) :*= s._1
      v(s._2, ::) :*= s._1
    })
    (u, v)
  }

}
