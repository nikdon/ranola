package ranola


import spire.implicits.cfor


trait RandomizedRangeFinder


/**
 * Algorithm 4.1 of "Finding structure with randomness:
 * Stochastic algorithms for constructing approximate matrix decompositions"
 * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
 */
object GenericRangeFinder extends RandomizedRangeFinder {

  /**
   * Find a projection matrix for the given one
   *
   * @param A           The input data matrix
   * @param sketchSize  The size of the matrix to return
   * @param op          The matrix operations (ex: multiplication etc.)
   * @return            A size-by-size projection matrix Q
   */
  def apply[N, M[_], V[_]](A: M[N], sketchSize: Int)(implicit op: MatrixOps[N, M, V]): M[N] = {
    val R = op.drawRandomMatrix(A, sketchSize)
    val Y = op.mltMM(A, op.mltMM(op.t(A), op.mltMM(A, R)))
    val (_q, _) = op.QR(Y)
    _q
  }
}


/**
 * Algorithm 4.3 of "Finding structure with randomness:
 * Stochastic algorithms for constructing approximate matrix decompositions"
 * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
 */
object PowerIterationsRangeFinder extends RandomizedRangeFinder {

  /**
   * Find a projection matrix for the given one
   *
   * @param A           The input data matrix
   * @param sketchSize  The size of the matrix to return
   * @param nIter       A number of power iterations used to stabilize the result
   * @param op          The matrix operations (ex: multiplication etc.)
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

    val (_q, _) = op.QR(Y)
    _q
  }
}


/**
 * Algorithm 4.4 of "Finding structure with randomness:
 * Stochastic algorithms for constructing approximate matrix decompositions"
 * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
 */
object SubspaceIterationsRangeFinder extends RandomizedRangeFinder {

  /**
   * Find a projection matrix for the given one
   *
   * @param A           The input data matrix
   * @param sketchSize  The size of the matrix to return
   * @param nIter       A number of power iterations used to stabilize the result
   * @param op          The matrix operations (ex: multiplication etc.)
   * @return            A size-by-size projection matrix Q
   */
  def apply[N, M[_], V[_]](A: M[N], sketchSize: Int, nIter: Int)(implicit op: MatrixOps[N, M, V]): M[N] = {
    var qi = op.drawZerosMatrix(op.getRows(A), math.min(math.min(op.getCols(A), op.getRows(A)), sketchSize))

    val R = op.drawRandomMatrix(A, sketchSize)
    val Y = op.mltMM(A, R)
    val (q, _) = op.QR(A)

    cfor(0)(_ < nIter, _ + 1){ i =>
      val Yih = op.mltMM(op.t(A), q)
      val (qih, _) = op.QR(Yih)
      val Yi = op.mltMM(A, qih)
      qi = op.QR(Yi)._1
    }

    qi
  }
}
