# ranola

**ranola** is a library for low-rank matrix approximations by randomized algorithms. Compared with standard deterministic algorithms, the randomized methods are often faster and produce accurate results.

**ranola** is built on top of [**scala breeze**][2] and utilizes big part of its resources.

At the moment the simplest direct scheme is implemented. It allows to ensure a minimal error in the final approximation of decompositions:

- Singular value decomposition:

      ```scala
      val SVD(u, s, v)
      ```

- Eigendecomposition:

      ```scala
      val EVD(l, v)
      ```

##  Randomized schemes

1. Generic scheme ([§4.1][1]) is designed for solving the fixed-rank problem, where the target rank of the input matrix is specified in advance. This algorithm works well for matrices whose singular values exhibit some decay, but they may produce a poor basis when the input matrix has a flat singular spectrum or when the input matrix is very large.

      ```scala
      val EVD(l, v) = evdr.generic(M: DenseMatrix[Double], k: Int, overSamples: Int)
      val SVD(u, s, v) = svdr.generic(M: DenseMatrix[Double], k: Int, overSamples: Int)
      ```
    
2. Adaptive scheme ([§4.2][1]) can handle the fixed-precision problem with predifined computational tolerance. The CPU time requirements of schemes 1 and 2 are essentially identical.

      ```scala
      val EVD(l, v) = evdr.adaptive(
        M: DenseMatrix[Double], 
        k: Int, tol: Double, 
        maxIter: Int, 
        overSamples: Int
      )
      
      val SVD(u, s, v) = svdr.adaptive(
        M: DenseMatrix[Double], 
        k: Int, 
        tol: Double, 
        maxIter: Int, 
        overSamples: Int
      )
      ```

3. Power Iteration scheme ([§4.3][1]) requires ```2*nIter + 1``` times as many matrix–vector multiplies as scheme 1, but is far more accurate in situations where the singular values of matrix to decompose decay slowly. Scheme 3 targets the fixed-rank problem. In situations where it is critical to achieve nearoptimal approximation errors, one can increase the oversampling beyond standard recommendation ```overSamples = 5``` all the way to ```overSamples = k``` without changing the scaling of the asymptotic computational cost. This procedure is vulnerable to round-off errors. The recommended implementation appears as scheme 4.

      ```scala
      val EVD(l, v) = evdr.powerIteration(
        M: DenseMatrix[Double], 
        k: Int, 
        nIter: Int, 
        overSamples: Int
      )
      
      val SVD(u, s, v) = svdr.powerIteration(
        M: DenseMatrix[Double], 
        k: Int, 
        nIter: Int, 
        overSamples: Int
      )
      ```

4. Subspace Iteration scheme ([§4.4][1]) is designed for solving the fixed-rank problem, where the target rank of the input matrix is specified in advance. This procedure is invulnerable to round-off errors in opposite to Power Iteration algorithm ([§4.3][1]).

      ```scala
      val EVD(l, v) = evdr.subspaceIteration(
        M: DenseMatrix[Double], 
        k: Int, 
        nIter: Int, 
        overSamples: Int
      )
        
      val SVD(u, s, v) = svdr.subspaceIteration(
        M: DenseMatrix[Double], 
        k: Int, 
        nIter: Int, 
        overSamples: Int
      )
      ```

5. Fast Generic scheme ([§4.5][1]) uses subsampled random Fourier transform (SRFT) that can be used to identify a nearoptimal basis for a rank-k matrix using l = (k + log(M.cols)) log(k) samples (```overSamples = (k + log(M.cols)) log(k) - k```). In practice, setting ```overSamples = 10``` or ```overSamples = 20``` is typically more than adequate. It allows to compute an approximate rank-(l) factorization of a general dense m-by-n matrix in roughly O(mn log(k + overSamples)) flops.

      ```scala
      val EVD(l, v) = evdr.fastGeneric(M: DenseMatrix[Double], k: Int, overSamples: Int)
      val SVD(u, s, v) = svdr.fastGeneric(M: DenseMatrix[Double], k: Int, overSamples: Int)
      ```

## Number of oversamples

Number of oversamples depends on:

- **Matrix dimensions**. Very large matrices may require more oversampling.

- **Singular spectrum**. The more rapid the decay of the singular values, the less oversampling is needed. In the extreme case that the matrix has exact rank k, it is not necessary to oversample.

- **Random test matrix**. Gaussian matrices succeed with very little oversampling, such as p = 5 or p = 10. There is rarely any advantage to select p > k.

## Notes

**ranola** is made based on ["Finding structure with randomness: Stochastic algorithms for constructing approximate matrix decompositions" Halko, et al., 2009 (arXiv:909)](http://arxiv.org/pdf/0909.4061)

[1]: http://arxiv.org/pdf/0909.4061
[2]: https://github.com/scalanlp/breeze
