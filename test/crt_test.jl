###############################################################################
###############################################################################
#
# This file contains units tests for the Chinese Remainder Theorem and it's 
# application to the multiplcation operation for PolynomialSparse128
#                                                                               
################################################################################
################################################################################

"""
Test Chinese Remainder Theorem implementation of multiplication.
"""
function crt_mult_test(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        @assert p1*p2 == chinese_mult(p1,p2)
    end
    println("crt_mult_test- PASSED")
end

#=need to do this properly 
"""
Test Chinese Remainder Theorem implementation of multiplication 
for polynomials with negative coefficients.
"""
function crt_mult_neg_test(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = -rand(PolynomialSparse128)
        p2 = -rand(PolynomialSparse128)
        println(p1)
        println(p2)
        @assert p1*p2 == chinese_mult(p1,p2)
    end
    println("crt_mult_test- PASSED")
end
=#