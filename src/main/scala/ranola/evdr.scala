package ranola


import breeze.linalg._
import breeze.linalg.eig.Eig
import breeze.linalg.eigSym.{DenseEigSym, EigSym}


/**
 * Approximate truncated randomized EVD
 */
object evdr {

  /**
   * @param M Matrix to decompose
   * @param k Number of eigenvalues and eigenvectors to extract
   * @param nOversamples Additional number of random vectors to sample the range of M so as
   *                     to ensure proper conditioning. The total number of random vectors
   *                     used to find the range of M is [k + nOversamples]
   * @return The eigenvalue decomposition (EVD) with the eigenvalues and the eigenvectors
   *
   *         ==References==
   *
   *         Finding structure with randomness: Stochastic algorithms for constructing
   *         approximate matrix decompositions
   *         Halko, et al., 2009 [[http://arxiv.org/abs/arXiv:0909.4061]]
   */
  def generic(M: DenseMatrix[Double], k: Int, nOversamples: Int): DenseEigSym = {
    require(k <= (M.rows min M.cols), "Number of columns in orthonormal matrix should be less than min(M.rows, M.cols)")
    require(k >= 1, "Sketch size should be greater than 1")

    val nRandom = k + nOversamples

    val Q = Generic(M, nRandom)

    val b = Q.t * (M * Q)

    val Eig(w, _, v) = eig(b)

    val _u = Q * v

    val u = utils.flipSigns(_u)

    EigSym(w(0 until k), u(::, 0 until k))
  }

  /**
   * @param M Matrix to decompose
   * @param k Number of eigenvalues and eigenvectors to extract
   * @param nOversamples Additional number of random vectors to sample the range of M so as
   *                     to ensure proper conditioning. The total number of random vectors
   *                     used to find the range of M is [k + nOversamples]
   * @param nIter Number of power iterations
   * @return The eigenvalue decomposition (EVD) with the eigenvalues and the eigenvectors
   *
   *         ==References==
   *
   *         Finding structure with randomness: Stochastic algorithms for constructing
   *         approximate matrix decompositions
   *         Halko, et al., 2009 [[http://arxiv.org/abs/arXiv:0909.4061]]
   */
  def powerIteration(M: DenseMatrix[Double], k: Int, nOversamples: Int, nIter: Int): DenseEigSym = {
    require(k <= (M.rows min M.cols), "Number of columns in orthonormal matrix should be less than min(M.rows, M.cols)")
    require(k >= 1, "Sketch size should be greater than 1")

    val nRandom = k + nOversamples

    val Q = PowerIteration(M, nRandom, nIter)

    val b = Q.t * (M * Q)

    val Eig(w, _, v) = eig(b)

    val _u = Q * v

    val u = utils.flipSigns(_u)

    EigSym(w(0 until k), u(::, 0 until k))
  }


}
