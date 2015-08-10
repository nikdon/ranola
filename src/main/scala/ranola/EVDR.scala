package ranola

import breeze.linalg.eig.Eig
import breeze.linalg.eigSym._
import breeze.linalg.{argmax, eig, DenseVector, DenseMatrix}
import breeze.numerics._


/**
 * Trait that describes eigenvalue decomposition in terms of randomized algorithms
 * Contains the only one function that resolves sign ambiguity
 *
 * @tparam N The type og values inside matrix and vector containers
 * @tparam M The matrix container
 * @tparam V The vector container
 */
trait EVDR[N, M[_], V[_], R] extends Decomposition[N, M, V, R] {

  /** Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive */
  def flipEVDSigns(u: M[N])(implicit op: MatrixOps[N, M, V]): M[N]
}


object EVDR {

  def viaPowerIteration[N, M[_], V[_], R](A: M[N], k: Int, nIter: Int, overSamples: Int)
                                         (implicit dec: EVDR[N, M, V, R],
                                          op: MatrixOps[N, M, V]) = {
    dec.viaPowerIteration(A, k, nIter, overSamples)
  }

  ///////////////////////////////////
  // Add new implementations below //
  ///////////////////////////////////

  implicit val breezeDoubleMatrixEVDR = new EVDR[Double, DenseMatrix, DenseVector, DenseEigSym] {

    override def decompose(A: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double])
                          (implicit op: MatrixOps[Double, DenseMatrix, DenseVector]): DenseEigSym = {
      val b = op.mltMM(op.t(Q), op.mltMM(A, Q))
      val Eig(w, _, v) = eig(b)
      val _u = op.mltMM(Q, v)
      val u = flipEVDSigns(_u)
      EigSym(w(0 until k), u(::, 0 until k))
    }

    override def flipEVDSigns(u: DenseMatrix[Double])
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