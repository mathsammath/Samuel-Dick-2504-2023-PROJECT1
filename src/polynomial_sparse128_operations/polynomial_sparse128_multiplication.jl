#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for PolynomialSparse128
#                                                                               
#############################################################################
#############################################################################
using Primes 

"""
Multiply two sparse polynomials.
"""
function *(p1::PolynomialSparse128, p2::PolynomialSparse128)::PolynomialSparse128
    p_out = PolynomialSparse128() #initialise result of multiplication
    p1_terms = collect(p1.lst)
    for t in p1_terms 
        p_out += t*p2 #multiply each term of p1 by p2
    end 
    return p_out 
end

"""
Power of a sparse polynomial.
"""
function ^(p::PolynomialSparse128, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

"""
Chinese Remainder Theorem Algorithm on two polynomials, 
both of type PolynomialSparse128.
Assume gcd(mod1, mod2) = 1
"""
#Based of pseudo-code presented in lectures.
#Should re-name variables.
function chinese_rem_thm(a::PolynomialSparse128, mod1::Int, b::PolynomialSparse128, mod2::Int)::PolynomialSparse128
    c = PolynomialSparse128()
    p1 = PolynomialModP(a, mod1)
    p2 = PolynomialModP(b, mod2)
    k = max(leading(p1).degree, leading(p2).degree)
    while k > -1
        if k ≤ leading(p1).degree
            ak = leading(p1).coeff # = coefficient on x^k in a     
            pop!(p1)
        else 
            ak = Int128(0) 
        end 
        if k ≤ leading(p2).degree
            bk = leading(p2).coeff # = coefficient on x^k in b
            pop!(p2)
        else 
            bk = Int128(0)
        end 
        ck = chinese_rem_thm([ak, bk], [mod1, mod2]) # integer chinese remainder
        c = c + ck*(x_poly_sparse128()^k)
        k -= 1
    end 
    return c
end 

"""
Multiplication for PolynomialSparse128 based on Chinese Remainder Theorem.
"""
function chinese_mult(a::PolynomialSparse128, b::PolynomialSparse128)::PolynomialSparse128
    height_a = findmax(coeffs(a))[1]
    height_b = findmax(coeffs(b))[1]

    B = 2*height_a*height_b*min(length(a), length(b))
    p = 3
    M = p
    c = mod(a*b, M)

    while M < B 
        p = nextprime(p, 2) #return next prime number 
        c_new = mod(a*b, p)
        c = chinese_rem_thm(c, M, c_new, p)
        M = M*p
    end 
    return mod(c, M)
end 

#=
def MULTIPLICATION: 
    Input   a = an x^n + ... a1 x + a0, b = bm x^m + ... b1 x + b0 in ZZ[x]
    Ouput   a * b in ZZ[x]

    height_a <- max(|an|, ..., |a1|, |a0|)
    height_b <- max(|bm|, ..., |b1|, |b0|)

    B <- 2 * height_a * height_b * min(n+1, m+1)
    p <- 3
    M <- p
    c <- (a*b) mod M          # multiplication in ZZ_3[x]

    while M < B do
       p  <- NextPrime(p)     # the first prime > p
       c' <- a*b mod p        # multiplication in ZZ_p[x]
       c  <- CRT([c, c'], [M, p])  # spatial optimaztion 
       M  <- M*p

    return c smod M           # restore negative coefficients with symemetric mod
=#

##########################################
# Helper functions
##########################################

"""
Chinese Remainder Theorem on Integers.
"""
function chinese_rem_thm(u::Vector{Int128}, m::Vector{Int})
    v = Vector{Int128}(undef, 2)
    v[1] = u[1]
    v[2] = (u[2] - v[1]) * Int128(int_inverse_mod(m[1], m[2]) % m[2])
    return v[1] + v[2]
end 