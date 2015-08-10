package ranola


/**
 * Trait that describes matrix decomposition in terms of randomized algorithms
 *
 * @tparam N The type og values inside matrix and vector containers
 * @tparam M The matrix container
 * @tparam V The vector container
 */
trait Decomposition[N, M[_], V[_], R] {

  /** Direct decomposition */
  protected [this] def decompose(A: M[N], k: Int, Q: M[N])(implicit op: MatrixOps[N, M, V]): R

  /**
   * Algorithm 4.3 of "Finding structure with randomness:
   * Stochastic algorithms for constructing approximate matrix decompositions"
   * Halko, et al., 2009 (arXiv:909) [[http://arxiv.org/pdf/0909.4061]]
   *
   * @param A     The input data matrix
   * @param k     The number of entries in result
   * @param nIter Number of iterations to stabilize the result
   * @return      A decomposition result
   */
  def viaPowerIteration(A: M[N], k: Int, nIter: Int, overSamples: Int)(implicit op: MatrixOps[N, M, V]): R = {
    val Q = PowerIterationRangeFinder(A, sketchSize = k + overSamples, nIter)
    decompose(A, k, Q)
  }
}



