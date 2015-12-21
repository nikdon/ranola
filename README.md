# ranola

[![Build Status](https://travis-ci.org/nikdon/ranola.svg?branch=master)](https://travis-ci.org/nikdon/ranola)
[![Codacy Badge](https://www.codacy.com/project/badge/134576e1957c4cd8a8cf1755bd839e71)](https://www.codacy.com/app/nd/ranola)
[ ![Download](https://api.bintray.com/packages/nikdon/ranola/ranola/images/download.svg) ](https://bintray.com/nikdon/ranola/ranola/_latestVersion)

**ranola** is a library for low-rank matrix approximations by randomized algorithms (Randomized Numerical Linear Algebra aka RandNLA). Compared with standard deterministic algorithms, the randomized methods are often faster and produce accurate results. Empirical results show that a corresponding implementation can be faster than LAPACK on dense matrices.

At the moment the simplest direct scheme is implemented. It allows to ensure a minimal error in the final approximation of decompositions:

- Singular value decomposition
- Eigen value decomposition

## Getting ranola

If you're using SBT, add the following line to your build file:

```scala
resolvers += "jitpack" at "https://jitpack.io"

libraryDependencies += "com.github.nikdon" % "ranola" % "v0.2.0"

```
    
For Maven:

```maven
<repository>
    <id>jitpack.io</id>
    <url>https://jitpack.io</url>
</repository>

<dependency>
    <groupId>com.github.nikdon</groupId>
    <artifactId>ranola</artifactId>
    <version>v0.2.0</version>
</dependency>

```

##  Randomized schemes

1. Generic scheme ([§4.1][1]) is designed for solving the fixed-rank problem, where the target rank of the input matrix is specified in advance. This algorithm works well for matrices whose singular values exhibit some decay, but they may produce a poor basis when the input matrix has a flat singular spectrum or when the input matrix is very large.

    ```scala
    val EigSym(lambda, evs) = EVDR.generic(A_dense, k = 3, nOverSamples = 2)
      
    val SVD(ur, sr, vr) = SVDR.generic(M, k = 3, nOverSamples = 10)
    
    ```
    
2. Adaptive scheme ([§4.2][1]) can handle the fixed-precision problem with predifined computational tolerance. The CPU time requirements of schemes 1 and 2 are essentially identical.

    ```scala
    not implemented in v 0.2.0
    ```

3. Power Iteration scheme ([§4.3][1]) requires ```2*nIter + 1``` times as many matrix–vector multiplies as scheme 1, but is far more accurate in situations where the singular values of matrix to decompose decay slowly. Scheme 3 targets the fixed-rank problem. In situations where it is critical to achieve nearoptimal approximation errors, one can increase the oversampling beyond standard recommendation ```overSamples = 5``` all the way to ```overSamples = k``` without changing the scaling of the asymptotic computational cost. This procedure is vulnerable to round-off errors. The recommended implementation appears as scheme 4.

    ```scala
    val SVD(ur, sr, vr) = SVDR.viaPowerIteration(M, k = 3, nIter = 5, nOverSamples = 2)
      
    val EigSym(lambda, evs) = EVDR.viaPowerIteration(M, k = 3, nIter = 5, nOverSamples = 2)
    
    ```

4. Subspace Iteration scheme ([§4.4][1]) is designed for solving the fixed-rank problem, where the target rank of the input matrix is specified in advance. This procedure is invulnerable to round-off errors in opposite to Power Iteration algorithm ([§4.3][1]).

    ```scala
    val SVD(ur, sr, vr) = SVDR.viaSubspaceIterations(M, k = 3, nIter = 5, nOverSamples = 2)
      
    val EigSym(lambda, evs) = EVDR.viaSubspaceIterations(M, k = 3, nIter = 5, nOverSamples = 2)
    
    ```

5. Fast Generic scheme ([§4.5][1]) uses subsampled random Fourier transform (SRFT) that can be used to identify a nearoptimal basis for a rank-k matrix using l = (k + log(M.cols)) log(k) samples (```overSamples = (k + log(M.cols)) log(k) - k```). In practice, setting ```overSamples = 10``` or ```overSamples = 20``` is typically more than adequate. It allows to compute an approximate rank-(l) factorization of a general dense m-by-n matrix in roughly O(mn log(k + overSamples)) flops.

    ```
    not implemented in v 0.2.0
    ```

## Number of oversamples

Number of oversamples depends on:

- **Matrix dimensions**. Very large matrices may require more oversampling.

- **Singular spectrum**. The more rapid the decay of the singular values, the less oversampling is needed. In the extreme case that the matrix has exact rank **k**, it is not necessary to oversample.

- **Random test matrix**. Gaussian matrices succeed with very little oversampling, such as ```overSamples = 5``` or ```overSamples = 10```. There is rarely any advantage to select **overSamples** > **k**.


## Add support of other math libraries

Anyone is welcome to customize **ranola**. The simplest way to it is add support of different math libraries. To do this find in ```MatrixOp.scala``` implicit values for custom matrix operations and add yours. For example, currently implemented the one one based on **scala breeze**:

```scala
implicit val BreezeDoubleMatrixOps = new MatrixOps[Double, breeze.linalg.DenseMatrix, breeze.linalg.DenseVector] {
    ...
}

```

Next step is to implement implicits in ```SVDR.scala``` and ```EVDR.scala```:

```scala
implicit val breezeDoubleMatrixSVDR = new Decomposition[Double, DenseMatrix, DenseVector, DenseSVD] {
    ...
}

implicit val breezeDoubleMatrixEVDR = new Decomposition[Double, DenseMatrix, DenseVector, DenseEigSym] {
...
}

```
    
To use the custom implementation simply import related implicits:

```scala
import ranola.EVDR.Implicits._

import ranola.SVDR.Implicits._

```

And that is it, enjoy.

[1]: http://arxiv.org/pdf/0909.4061
