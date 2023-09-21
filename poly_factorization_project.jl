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
include("src/term128.jl")
include("src/polynomial_dense.jl")
    include("src/polynomial_dense_operations/polynomial_dense_addition.jl")
    include("src/polynomial_dense_operations/polynomial_dense_multiplication.jl")
    include("src/polynomial_dense_operations/polynomial_dense_division.jl")
    include("src/polynomial_dense_operations/polynomial_dense_gcd.jl")
    include("src/polynomial_factorization/factor.jl")

include("src/polynomial_sparse.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_addition.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_multiplication.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_division.jl")
    include("src/polynomial_sparse_operations/polynomial_sparse_gcd.jl")

include("src/polynomial_sparse128.jl")
    include("src/polynomial_sparse128_operations/polynomial_sparse128_addition.jl")
    include("src/polynomial_sparse128_operations/polynomial_sparse128_multiplication.jl")
    include("src/polynomial_sparse128_operations/polynomial_sparse128_division.jl")
    include("src/polynomial_sparse128_operations/polynomial_sparse128_gcd.jl")

include("src/polynomial_modp.jl")
    include("src/polynomial_modp_operations/polynomial_modp_addition.jl")
    include("src/polynomial_modp_operations/polynomial_modp_division.jl")
    include("src/polynomial_modp_operations/polynomial_modp_gcd.jl")
    include("src/polynomial_modp_operations/polynomial_modp_multiplication.jl")

include("test/crt_test.jl")
include("test/factorization_test.jl")
include("test/integers_test.jl")
include("test/polynomials_test.jl")
include("test/polynomial_modp_test.jl")
include("test/polynomial_sparse128_test.jl")
include("test/polynomial_sparse_test.jl")
include("test/pow_mod.test.jl")
nothing