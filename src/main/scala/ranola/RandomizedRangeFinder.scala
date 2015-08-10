package ranola


trait RandomizedRangeFinder


/**
 * Algorithm 4.3 of "Finding structure with randomness:
 * Stochastic algorithms for constructing approximate matrix decompositions"
 * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
 */
object PowerIterationRangeFinder extends RandomizedRangeFinder {

  /**
   * Find a projection matrix for the given one
   *
   * @param A           The input data matrix
   * @param sketchSize  Size of the matrix to return
   * @param nIter       Number of power iterations used to stabilize the result
   * @return            A size-by-size projection matrix Q
   */
  def apply[N, M[_], V[_]](A: M[N], sketchSize: Int, nIter: Int)(implicit op: MatrixOps[N, M, V]): M[N] = {
    val R = op.drawRandomMatrix(A, sketchSize)
    var Y = op.mltMM(A, R)

    var i = 0
    while (i < nIter) {
      Y = op.mltMM(A, op.mltMM(op.t(A), Y))
      i += 1
    }

    val (_q, _) = op.QR(A)
    _q
  }
}
