package ranola


import breeze.linalg.svd.{DenseSVD, SVD}
import breeze.linalg.{DenseMatrix, min, svd}


/**
 * Computes an approximate truncated randomized SVD. Fast on large matrices with limited
 * number of singular values and vectors that are necessary to extract
 */
object svdr {


  object generic {
    def apply(M: DenseMatrix[Double], k: Int, overSamples: Int): DenseSVD = {
      doSvdr(M, k, RandomizedRangeFinder.generic(M, sketchSize = k + overSamples))
    }
  }


  object fastGeneric {
    def apply(M: DenseMatrix[Double], k: Int, overSamples: Int): DenseSVD = {
      doSvdr(M, k, RandomizedRangeFinder.fastGeneric(M, sketchSize = k + overSamples))
    }
  }


  object powerIteration {
    def apply(M: DenseMatrix[Double], k: Int, nIter: Int, overSamples: Int): DenseSVD = {
      doSvdr(M, k, RandomizedRangeFinder.powerIteration(M, sketchSize = k + overSamples, nIter))
    }
  }


  object adaptive {
    def apply(M: DenseMatrix[Double], k: Int, tol: Double, maxIter: Int, overSamples: Int): DenseSVD = {
      doSvdr(M, k, RandomizedRangeFinder.adaptive(M, nRandVec = k + overSamples, tol, maxIter))
    }
  }


  object subspaceIteration {
    def apply(M: DenseMatrix[Double], k: Int, nIter: Int, overSamples: Int): DenseSVD = {
      doSvdr(M, k, RandomizedRangeFinder.subspaceIteration(M, sketchSize = k + overSamples, nIter))
    }
  }


  private def doSvdr(M: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double]): DenseSVD = {
    require(k <= min(Q.rows, Q.cols), "min(Q.rows, Q.cols) should be less or equal to k")

    val b = Q.t * M
    val SVD(w2, _s, _v) = svd.reduced(b)
    val _u = Q * w2
    val (u, v) = utils.flipSVDSigns(_u, _v)
    SVD(u(::, 0 until k), _s(0 until k), v(0 until k, ::))
  }
}
