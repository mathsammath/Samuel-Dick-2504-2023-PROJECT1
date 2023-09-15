#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Distributions, StatsBase, Random

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

include("src/general_alg.jl")
include("src/term.jl")
include("src/polynomial_dense.jl")
    include("src/polynomial_dense_operations/polynomial_dense_addition.jl")
    include("src/polynomial_dense_operations/polynomial_dense_multiplication.jl")
    include("src/polynomial_dense_operations/polynomial_dense_division.jl")
    include("src/polynomial_dense_operations/polynomial_dense_gcd.jl")
include("src/polynomial_dense_factorization/factor.jl")

include("src/polynomial_sparse.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_addition.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_multiplication.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_division.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_gcd.jl")

nothing