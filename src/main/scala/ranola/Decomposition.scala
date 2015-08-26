package ranola


/**
 * Trait that describes matrix decomposition in terms of randomized algorithms
 *
 * @tparam N The type of values inside matrix and vector containers
 * @tparam M The matrix container
 * @tparam V The vector container
 * @tparam R The result container
 */
trait Decomposition[N, M[_], V[_], R] {

  /**
   * Direct decomposition
   *
   * @param A   The matrix to decompose
   * @param k   The number of entries in the result
   * @param Q   Orthonormal matrix of M
   * @param op  Matrix operations
   * @return    A result of decomposition
   */
  protected [this] def decompose(A: M[N], k: Int, Q: M[N])(implicit op: MatrixOps[N, M, V]): R

  /**
   * Algorithm 4.1 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param A             The input data matrix
   * @param k             The number of entries in result
   * @param nOverSamples  The number of oversamples
   * @param op            The matrix operations (ex: multiplication etc.)
   * @return              A decomposition result
   */
  def generic(A: M[N], k: Int, nOverSamples: Int)(implicit op: MatrixOps[N, M, V]): R = {
    val Q = GenericRangeFinder(A, sketchSize = k + nOverSamples)
    decompose(A, k, Q)
  }

  /**
   * Algorithm 4.3 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param A             The input data matrix
   * @param k             The number of entries in result
   * @param nIter         Number of iterations to stabilize the result
   * @param nOverSamples  The number of oversamples
   * @param op            The matrix operations (ex: multiplication etc.)
   * @return              A decomposition result
   */
  def viaPowerIterations(A: M[N], k: Int, nIter: Int, nOverSamples: Int)(implicit op: MatrixOps[N, M, V]): R = {
    val Q = PowerIterationsRangeFinder(A, sketchSize = k + nOverSamples, nIter)
    decompose(A, k, Q)
  }

  /**
   * Algorithm 4.4 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param A             The input data matrix
   * @param k             The number of entries in result
   * @param nIter         Number of iterations to stabilize the result
   * @param nOverSamples  The number of oversamples
   * @param op            The matrix operations (ex: multiplication etc.)
   * @return              A decomposition result
   */
  def viaSubspaceIterations(A: M[N], k: Int, nIter: Int, nOverSamples: Int)(implicit op: MatrixOps[N, M, V]): R = {
    val Q = SubspaceIterationsRangeFinder(A, sketchSize = k + nOverSamples, nIter)
    decompose(A, k, Q)
  }
}



