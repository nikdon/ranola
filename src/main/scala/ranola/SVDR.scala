package ranola

import breeze.linalg._
import breeze.linalg.svd.{DenseSVD, SVD}
import breeze.numerics._


object SVDR {

  def generic[N, M[_], V[_], R](A: M[N], k: Int, nOverSamples: Int)
                               (implicit dec: Decomposition[N, M, V, R], op: MatrixOps[N, M, V]) = {
    dec.generic(A, k, nOverSamples)
  }

  def viaPowerIteration[N, M[_], V[_], R](A: M[N], k: Int, nIter: Int, nOverSamples: Int)
                                         (implicit dec: Decomposition[N, M, V, R], op: MatrixOps[N, M, V]) = {
    dec.viaPowerIteration(A, k, nIter, nOverSamples)
  }


  object Implicits {

    //////////////////////////////////
    // Add new implementations here //
    //////////////////////////////////

    implicit val breezeDoubleMatrixSVDR = new Decomposition[Double, DenseMatrix, DenseVector, DenseSVD] {

      override def decompose(A: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double])
                            (implicit op: MatrixOps[Double, DenseMatrix, DenseVector]): DenseSVD = {
        val b = op.mltMM(op.t(Q), A)
        val SVD(w2, _s, _v) = svd.reduced(b)
        val _u = op.mltMM(Q, w2)
        val (u, v) = flipSVDSigns(_u, _v)
        SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
      }

      def flipSVDSigns(u: DenseMatrix[Double], v: DenseMatrix[Double])
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

}