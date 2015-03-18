package ranola


import breeze.linalg.svd.{DenseSVD, SVD}
import breeze.linalg.{DenseMatrix, svd}


/**
 * Computes an approximate truncated randomized SVD. Fast on large matrices with limited
 * number of singular values and vectors that are necessary to extract
 */
object svdr {

  def apply(M: DenseMatrix[Double], Q: DenseMatrix[Double], k: Int): DenseSVD = {
    doSvdr(M, Q, k)
  }

  private def doSvdr(M: DenseMatrix[Double], Q: DenseMatrix[Double], k: Int): DenseSVD = {
    require(k <= Q.cols, "Number of columns Q should be less or equal to k")

    val b = Q.t * M
    val SVD(w2, _s, _v) = svd.reduced(b)
    val _u = Q * w2
    val (u, v) = utils.flipSVDSigns(_u, _v)
    SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
  }
}
