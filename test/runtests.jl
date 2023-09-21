#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")
####
# Execute unit tests for integers
###
include("integers_test.jl")
test_euclid_ints()
test_ext_euclid_ints()

####
# Execute unit tests for polynomials (type PolynomialDense)
####
include("polynomials_test.jl")
prod_test_poly()
prod_derivative_test_poly()
ext_euclid_test_poly()
division_test_poly()

####
# Execute unit tests for polynomials (type PolynomialSparse)
####
include("polynomial_sparse_test.jl")
prod_test_poly_sparse()
prod_derivative_test_poly_sparse()
division_test_poly_sparse()
ext_euclid_test_poly_sparse()
push_test_poly_sparse()
pop_test_poly_sparse()
mod_test_poly_sparse()

####
# Execute unit tests for polynomials (type PolynomialSparse128)
####
include("polynomial_sparse128_test.jl")
prod_test_poly_sparse128()
prod_derivative_test_poly_sparse128()
division_test_poly_sparse128()
ext_euclid_test_poly_sparse128()
push_test_poly_sparse128()
pop_test_poly_sparse128()
mod_test_poly_sparse128()
bigint_test()

####
# Execute unit tests for polynomials (type PolynomialModP)
####
include("polynomial_modp_test.jl")
prod_test_poly_modp()
prod_derivative_test_poly_modp()
div_rem_test_poly_modp()
push_test_poly_modp()
pop_test_poly_modp()


####
# Execute unit tests for Chinese Remainder Theorem and associated multiplication
####
include("crt_test.jl")
crt_mult_test()

####
# Execute unit tests for ^ and pow_mod refactored methods.
####
include("pow_mod.test.jl")
pow_sparse_test()
pow_sparse128_test()
powmod_sparse_test()
powmod_sparse128_test()


####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
factor_test_poly()

