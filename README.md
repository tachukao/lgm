# LGM: lightweight package for factor analysis and probabilistic PCA

This package is built on top of [Owl](https://github.com/owlbarn/owl.git). The EM algorithm is used to do inference for both factor analysis (FA) and probabilistic PCA (PPCA).

**Install:**
```sh
make install
```

**Uninstall:**
```sh
make uninstall
```

**Usage:**

FA
```ocaml
#require "lgm"
let mu, w, sigma, zs, nll = Lgm.infer ~model:`fa data d 
```

PPCA
```ocaml
#require "lgm"
let mu, w, sigma, zs, nll = Lgm.infer ~model:`ppca data d 
```

Here, `mu` is the mean of the data set (dims `n x 1`).
`w` is the loading matrix (dims `n x d`). 
`sigma` is the noise covariance (dims `n x n`). 
`zs` is a matrix of inferred latent variables (dims `d x n_samples`). 
`data` is a data matrix (dims `n x n_samples`). 
`nll` is the negative log likelihood of data for the inferred parameters. 
