package ranola

import breeze.linalg._
import breeze.linalg.eig.Eig
import breeze.linalg.eigSym.{DenseEigSym, EigSym}
import breeze.linalg.support.CanTranspose
import breeze.linalg.svd.{DenseSVD, SVD}
import breeze.numerics._


trait Decomposition {

  type Result

  /**
   * Abstract decomposition method
   *
   * @param M Matrix to decompose
   * @param k Number of entries in result
   * @param Q Orthonormal matrix of M
   * @return Result of decomposition
   */
  def decompose[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, Q: DenseMatrix[Double])
                                                (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                 mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                 trans: CanTranspose[M, MT],
                                                 mltTrans: OpMulMatrixDenseMatrix[MT])
  : Result

  def generic[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, overSamples: Int)
                                              (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                               mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                               trans: CanTranspose[M, MT],
                                               mltTrans: OpMulMatrixDenseMatrix[MT]) = {
    decompose(M, k, GenericRangeFinder(M, sketchSize = k + overSamples))
  }

  def fastGeneric[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, overSamples: Int)
                                                  (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                   mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                   trans: CanTranspose[M, MT],
                                                   mltTrans: OpMulMatrixDenseMatrix[MT]) = {
    decompose(M, k, FastGenericRangeFinder(M, sketchSize = k + overSamples))
  }

  def powerIteration[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, nIter: Int, overSamples: Int)
                                                     (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                      mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                      trans: CanTranspose[M, MT],
                                                      mltTrans: OpMulMatrixDenseMatrix[MT]) = {
    decompose(M, k, PowerIterationRangeFinder(M, sketchSize = k + overSamples, nIter))
  }

  def subspaceIteration[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, nIter: Int, overSamples: Int)
                                                        (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                         mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                         trans: CanTranspose[M, MT],
                                                         mltTrans: OpMulMatrixDenseMatrix[MT]) = {
    decompose(M, k, SubspaceIterationRangeFinder(M, sketchSize = k + overSamples, nIter))
  }

  def adaptive[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, tol: Double, maxIter: Int, overSamples: Int)
                                               (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                trans: CanTranspose[M, MT],
                                                mltTrans: OpMulMatrixDenseMatrix[MT],
                                                mltMatDenVec: OpMulMatrixDenseVector[M]) = {
    decompose(M, k, AdaptiveRangeFinder(M, nRandVec = k + overSamples, tol, maxIter))
  }
}


object evdr extends Decomposition {

  type Result = DenseEigSym

  /** Direct randomized eigendecomposition */
  override def decompose[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, Q: DenseMatrix[Double])
                                                         (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                          mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                          trans: CanTranspose[M, MT],
                                                          mltTrans: OpMulMatrixDenseMatrix[MT])
  : Result = {

    require(k <= min(Q.rows, Q.cols), "min(Q.rows, Q.cols) should be less or equal to k")

    val b = Q.t * mltMatDenMat(M, Q) // Q.t * (M * Q)
    val Eig(w, _, v) = eig(b)
    val _u = Q * v
    val u = flipEVDSigns(_u)
    EigSym(w(0 until k), u(::, 0 until k))
  }

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive.
   *
   * @param u eigenvectors
   * @return eigenvectors with resolved sign ambiguity
   */
  private[this] def flipEVDSigns(u: DenseMatrix[Double]): DenseMatrix[Double] = {
    val abs_u = abs(u)
    val len = u.cols
    var i = 0
    while (i < len) {
      u(::, i) :*= signum(u(argmax(abs_u(::, i)), i))
      i += 1
    }
    u
  }
}


object svdr extends Decomposition {

  type Result = DenseSVD

  /** Direct randomized singular value decomposition */
  override def decompose[M <: AnyMatrix, MT <: AnyMatrix](M: M, k: Int, Q: DenseMatrix[Double])
                                                         (implicit mltMatDenMat: OpMulMatrixDenseMatrix[M],
                                                          mltDenMatMat: OpMulDenseMatrixMatrix[M],
                                                          trans: CanTranspose[M, MT],
                                                          mltTrans: OpMulMatrixDenseMatrix[MT])
  : Result = {

    require(k <= min(Q.rows, Q.cols), "min(Q.rows, Q.cols) should be less or equal to k")

    val b = mltDenMatMat(Q.t, M) // Q.t * M
    val SVD(w2, _s, _v) = svd.reduced(b)
    val _u = Q * w2
    val (u, v) = flipSVDSigns(_u, _v)
    SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
  }

  /**
   * Resolves the sign ambiguity. Largest in absolute value entries of u columns are always positive.
   *
   * @param u left singular vectors
   * @param v right singular vectors
   * @return left and right singular vectors with resolved sign ambiguity
   */
  private[this] def flipSVDSigns(u: DenseMatrix[Double],
                                 v: DenseMatrix[Double]): (DenseMatrix[Double], DenseMatrix[Double]) = {
    val abs_u = abs(u)
    val len = u.cols
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
