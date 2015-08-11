package ranola

import breeze.linalg._
import breeze.linalg.svd.{DenseSVD, SVD}
import breeze.numerics._


/**
 * Trait that describes singular value decomposition in terms of randomized algorithms
 * Contains the only one function that resolves sign ambiguity
 *
 * @tparam N The type of values inside matrix and vector containers
 * @tparam M The matrix container
 * @tparam V The vector container
 * @tparam R The result container
 */
trait SVDR[N, M[_], V[_], R] extends Decomposition[N, M, V, R] {

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive
   *
   * @param u   Left singular vectors
   * @param v   Right singular vectors
   * @param op  Matrix operations
   * @return    Left and right singular vectors with resolved sign ambiguity
   */
  def flipSVDSigns(u: M[N], v: M[N])(implicit op: MatrixOps[N, M, V]): (M[N], M[N])
}


object SVDR {

  def generic[N, M[_], V[_], R](A: M[N], k: Int, nOverSamples: Int)
                               (implicit dec: SVDR[N, M, V, R],
                                op: MatrixOps[N, M, V]) = {
    dec.generic(A, k, nOverSamples)
  }

  def viaPowerIteration[N, M[_], V[_], R](A: M[N], k: Int, nIter: Int, nOverSamples: Int)
                                         (implicit dec: SVDR[N, M, V, R],
                                          op: MatrixOps[N, M, V]) = {
    dec.viaPowerIteration(A, k, nIter, nOverSamples)
  }

  ///////////////////////////////////
  // Add new implementations below //
  ///////////////////////////////////

  implicit val breezeDoubleMatrixSVDR = new SVDR[Double, DenseMatrix, DenseVector, DenseSVD] {

    override def decompose(A: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double])
                          (implicit op: MatrixOps[Double, DenseMatrix, DenseVector]): DenseSVD = {
      val b = op.mltMM(op.t(Q), A)
      val SVD(w2, _s, _v) = svd.reduced(b)
      val _u = op.mltMM(Q, w2)
      val (u, v) = flipSVDSigns(_u, _v)
      SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
    }

    override def flipSVDSigns(u: DenseMatrix[Double], v: DenseMatrix[Double])
                             (implicit op: MatrixOps[Double, DenseMatrix, DenseVector]): (DenseMatrix[Double], DenseMatrix[Double]) = {
      val abs_u = op.absM(u)
      val len = op.getCols(u)
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

}