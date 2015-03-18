package ranola


import breeze.linalg.eig.Eig
import breeze.linalg.eigSym.{DenseEigSym, EigSym}
import breeze.linalg.{DenseMatrix, eig}


/**
 * Computes an approximate truncated randomized EVD. Fast on large matrices with limited
 * number of singular values and vectors that are necessary to extract
 */
object evdr {

  def apply(M: DenseMatrix[Double], Q: DenseMatrix[Double], k: Int): DenseEigSym = {
    val b = Q.t * (M * Q)
    val Eig(w, _, v) = eig(b)
    val _u = Q * v
    val u = utils.flipSigns(_u)
    EigSym(w(0 until k), u(::, 0 until k))
  }
}
