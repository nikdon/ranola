package ranola

import breeze.linalg.eig.Eig
import breeze.linalg.eigSym._
import breeze.linalg.{DenseMatrix, DenseVector, argmax, eig}
import breeze.numerics._


object EVDR {

  def generic[N, M[_], V[_], R](A: M[N], k: Int, nOverSamples: Int)
                               (implicit dec: Decomposition[N, M, V, R], op: MatrixOps[N, M, V]) = {
    dec.generic(A, k, nOverSamples)
  }

  def viaPowerIterations[N, M[_], V[_], R](A: M[N], k: Int, nIter: Int, nOverSamples: Int)
                                         (implicit dec: Decomposition[N, M, V, R], op: MatrixOps[N, M, V]) = {
    dec.viaPowerIterations(A, k, nIter, nOverSamples)
  }

  def viaSubspaceIterations[N, M[_], V[_], R](A: M[N], k: Int, nIter: Int, nOverSamples: Int)
                                         (implicit dec: Decomposition[N, M, V, R], op: MatrixOps[N, M, V]) = {
    dec.viaSubspaceIterations(A, k, nIter, nOverSamples)
  }


  object Implicits {

    //////////////////////////////////
    // Add new implementations here //
    //////////////////////////////////

    implicit val breezeDoubleMatrixEVDR = new Decomposition[Double, DenseMatrix, DenseVector, DenseEigSym] {

      override def decompose(A: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double])
                            (implicit op: MatrixOps[Double, DenseMatrix, DenseVector]): DenseEigSym = {
        val b = op.mltMM(op.t(Q), op.mltMM(A, Q))
        val Eig(w, _, v) = eig(b)
        val _u = op.mltMM(Q, v)
        val u = flipEVDSigns(_u)
        EigSym(w(0 until k), u(::, 0 until k))
      }

      def flipEVDSigns(u: DenseMatrix[Double])
                      (implicit op: MatrixOps[Double, DenseMatrix, DenseVector]): DenseMatrix[Double] = {
        val abs_u = op.absM(u)
        val len = op.getCols(u)
        var i = 0
        while (i < len) {
          u(::, i) :*= signum(u(argmax(abs_u(::, i)), i))
          i += 1
        }
        u
      }
    }
  }

}