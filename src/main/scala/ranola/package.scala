import breeze.linalg.operators.OpMulMatrix
import breeze.linalg.{DenseMatrix, DenseVector, Matrix}

package object ranola {

  type AnyMatrix = Matrix[_]

  // Multiplications
  type OpMulMatrixDenseMatrix[M] = OpMulMatrix.Impl2[M, DenseMatrix[Double], DenseMatrix[Double]]
  type OpMulMatrixDenseVector[M] = OpMulMatrix.Impl2[M, DenseVector[Double], DenseVector[Double]]
  type OpMulDenseMatrixMatrix[M] = OpMulMatrix.Impl2[DenseMatrix[Double], M, DenseMatrix[Double]]

}
