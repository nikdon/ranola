package ranola

import breeze.linalg._
import breeze.linalg.eig.Eig
import breeze.linalg.eigSym.{DenseEigSym, EigSym}
import breeze.linalg.svd.{DenseSVD, SVD}


trait Decomposition[T] {

  def decompose(M: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double]): T

  def generic(M: DenseMatrix[Double], k: Int, overSamples: Int) = {
    decompose(M, k, GenericRangeFinder(M, sketchSize = k + overSamples))
  }

  def fastGeneric(M: DenseMatrix[Double], k: Int, overSamples: Int) = {
    decompose(M, k, FastGenericRangeFinder(M, sketchSize = k + overSamples))
  }

  def powerIteration(M: DenseMatrix[Double], k: Int, nIter: Int, overSamples: Int) = {
    decompose(M, k, PowerIterationRangeFinder(M, sketchSize = k + overSamples, nIter))
  }

  def adaptive(M: DenseMatrix[Double], k: Int, tol: Double, maxIter: Int, overSamples: Int) = {
    decompose(M, k, AdaptiveRangeFinder(M, nRandVec = k + overSamples, tol, maxIter))
  }

  def subspaceIteration(M: DenseMatrix[Double], k: Int, nIter: Int, overSamples: Int) = {
    decompose(M, k, SubspaceIterationRangeFinder(M, sketchSize = k + overSamples, nIter))
  }
}


object evdr extends Decomposition[DenseEigSym] {

  override def decompose(M: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double]): DenseEigSym = {
    require(k <= min(Q.rows, Q.cols), "min(Q.rows, Q.cols) should be less or equal to k")

    val b = Q.t * (M * Q)
    val Eig(w, _, v) = eig(b)
    val _u = Q * v
    val u = utils.flipEVDSigns(_u)
    EigSym(w(0 until k), u(::, 0 until k))
  }
}


object svdr extends Decomposition[DenseSVD] {

  override def decompose(M: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double]): DenseSVD = {
    require(k <= min(Q.rows, Q.cols), "min(Q.rows, Q.cols) should be less or equal to k")

    val b = Q.t * M
    val SVD(w2, _s, _v) = svd.reduced(b)
    val _u = Q * w2
    val (u, v) = utils.flipSVDSigns(_u, _v)
    SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
  }
}
