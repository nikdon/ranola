package ranola


import breeze.linalg.eig.Eig
import breeze.linalg.eigSym.{DenseEigSym, EigSym}
import breeze.linalg.{DenseMatrix, eig, min}


/**
 * Computes an approximate truncated randomized EVD. Fast on large matrices with limited
 * number of singular values and vectors that are necessary to extract
 */
object evdr {


  object generic {
    def apply(M: DenseMatrix[Double], k: Int, overSamples: Int): DenseEigSym = {
      doEvdr(M, k, RandomizedRangeFinder.generic(M, sketchSize = k + overSamples))
    }
  }


  object fastGeneric {
    def apply(M: DenseMatrix[Double], k: Int, overSamples: Int): DenseEigSym = {
      doEvdr(M, k, RandomizedRangeFinder.fastGeneric(M, sketchSize = k + overSamples))
    }
  }


  object powerIteration {
    def apply(M: DenseMatrix[Double], k: Int, nIter: Int, overSamples: Int): DenseEigSym = {
      doEvdr(M, k, RandomizedRangeFinder.powerIteration(M, sketchSize = k + overSamples, nIter))
    }
  }


  object adaptive {
    def apply(M: DenseMatrix[Double], k: Int, tol: Double, maxIter: Int, overSamples: Int): DenseEigSym = {
      doEvdr(M, k, RandomizedRangeFinder.adaptive(M, nRandVec = k + overSamples, tol, maxIter))
    }
  }


  object subspaceIteration {
    def apply(M: DenseMatrix[Double], k: Int, nIter: Int, overSamples: Int): DenseEigSym = {
      doEvdr(M, k, RandomizedRangeFinder.subspaceIteration(M, sketchSize = k + overSamples, nIter))
    }
  }


  private def doEvdr(M: DenseMatrix[Double], k: Int, Q: DenseMatrix[Double]): DenseEigSym = {
    require(k <= min(Q.rows, Q.cols), "min(Q.rows, Q.cols) should be less or equal to k")

    val b = Q.t * (M * Q)
    val Eig(w, _, v) = eig(b)
    val _u = Q * v
    val u = utils.flipEVDSigns(_u)
    EigSym(w(0 until k), u(::, 0 until k))
  }
}
