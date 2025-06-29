# 2504_2023_project1 - Polynomial Factorization

This project implements polynomial arithmetic and polynomial factorization for polynomials with integer coefficients. 

To load all functionality, in the directory of the package:

```
] activate .
```

```
julia> include("poly_factorization_project.jl")
```

You may then use functionality such as,

```
julia> gcd(rand(Polynomial) + rand(Polynomial), rand(Polynomial), 101)
```

To execute all unit tests run:

```
julia> include("test/runtests.jl")
```

You may see examples in `example_script.jl` and run that script line by line.
