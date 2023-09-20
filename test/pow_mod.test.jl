###############################################################################
###############################################################################
#
# This file contains units tests for the new ^ and pow_mod methods for  
# polynomials of type PolynomialSparse and PolynomialSparse128.
#                                                                               
################################################################################
################################################################################

using Primes 

"""
^ method using repeated squaring test for PolynomialSparse types.
"""
function pow_sparse_test(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        n = rand(1:5) #for larger n, Sparse type overflows
        #println(p1)
        #println(n)
        @assert p1^n == old_pow(p1, n)
    end
    println("pow_sparse_test- PASSED")
end

"""
^ method using repeated squaring test for PolynomialSparse128 types.
"""
function pow_sparse128_test(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        n = rand(1:5)
        @assert p1^n == old_pow(p1, n)
    end
    println("pow_sparse128_test- PASSED")
end

"""
pow_mod method using repeated squaring test for PolynomialSparse types.
"""
function powmod_sparse_test(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        n = rand(1:5) #for larger n, Sparse type overflows
        p = rand(primes(1,20))
        @assert pow_mod(p1, n, p) == old_pow_mod(p1, n, p)
    end
    println("powmod_sparse_test- PASSED")
end

"""
pow_mod method using repeated squaring test for PolynomialSparse128 types.
"""
function powmod_sparse128_test(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        n = rand(1:20)
        p = rand(primes(1,20))
        @assert pow_mod(p1, n, p) == old_pow_mod(p1, n, p)
    end
    println("crt_mult_test- PASSED")
end




##################################################
# Old implementations of ^ and pow_mod methods
##################################################
"""
Power of a sparse polynomial.
"""
function old_pow(p::PolynomialSparse128, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end 

"""
Power of a sparse polynomial.
"""
function old_pow(p::PolynomialSparse, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end 

"""
Power of a sparse polynomial mod prime. Re-factored using repeated squaring.
"""
function old_pow_mod(p::PolynomialSparse, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end 

"""
Power of a sparse polynomial mod prime.
"""
function old_pow_mod(p::PolynomialSparse128, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end 