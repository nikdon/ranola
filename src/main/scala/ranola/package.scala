package object ranola {

  import breeze.linalg.operators.OpMulMatrix
  import breeze.linalg.{DenseMatrix, DenseVector, Matrix}

  type AnyMatrix = Matrix[_]

  // Multiplications
  type OpMulMatrixDenseMatrix[M <: AnyMatrix] = OpMulMatrix.Impl2[M, DenseMatrix[Double], DenseMatrix[Double]]
  type OpMulMatrixDenseVector[M <: AnyMatrix] = OpMulMatrix.Impl2[M, DenseVector[Double], DenseVector[Double]]
  type OpMulDenseMatrixMatrix[M <: AnyMatrix] = OpMulMatrix.Impl2[DenseMatrix[Double], M, DenseMatrix[Double]]

}
