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
Power of a polynomial (type Sparse128).
"""
function ^(p::PolynomialSparse128, n::Int)
    n < 0 && error("No negative power")
    b = reverse(string(n; base=2)) #binary string representation of n, reversed.
    ans, w = Int128(1), p
    for i in 1:length(b)
        if b[i] == '1'
            ans = ans*w
        end
        w = w*w
    end 
    return ans 
end

"""
Chinese Remainder Theorem Algorithm on two polynomials, 
both of type PolynomialSparse128. Assume gcd(mod1, mod2) = 1.
"""
function chinese_rem_thm(p::Vector{PolynomialSparse128}, mods::Vector{Int})::PolynomialSparse128
    p1 = PolynomialModP(p[1], mods[1]) #instantiate polynomials as type PolynomialModP
    p2 = PolynomialModP(p[2], mods[2])
    result = zero(PolynomialSparse128) #initialise result of algorithm
    k = max(degree(p1), degree(p2)) 
    while k > -1 
        if k > degree(p1) #no xᵏ term exists in p1
            ak = Int128(0)
        else 
            ak = Int128(leading(p1).coeff) #coefficient of xᵏ of p1 
            pop!(p1)
        end 
        if k > degree(p2) #no xᵏ term exists in p2
             bk = Int128(0) 
        else 
            bk = Int128(leading(p2).coeff) #coefficinet of xᵏ of p2
            pop!(p2)
        end 
        ck = chinese_rem_thm([ak, bk], [p1.prime_mod, p2.prime_mod])
        result = result + ck*(x_poly_sparse128()^k)
        k = k - 1
    end 
    return result 
end 

"""
Multiplication for PolynomialSparse128 based on Chinese Remainder Theorem.
Code adapted from algorithm provided in lectures.
"""
function chinese_mult(a::PolynomialSparse128, b::PolynomialSparse128)::PolynomialSparse128
    height_a = findmax(coeffs(a))[1]
    height_b = findmax(coeffs(b))[1]

    B = 2*height_a*height_b*min(length(a)+1, length(b)+1) #bound 
    p = 3 #first prime 
    M = p 
    result = mod(a*b, M) #initialise eventual result of mult

    while M < B 
        p = nextprime(p, 2) #return next prime number 
        c_new = mod(a*b, p)
        result = chinese_rem_thm([result, c_new], [M, p])
        M = M*p
    end 
    return mod(result, M)
end 

##########################################
# Helper functions
##########################################

"""
Chinese Remainder Theorem on two Integers.
"""
function chinese_rem_thm(u::Vector{Int128}, m::Vector{Int})
    v = Vector{Int128}(undef, 2)
    v[1] = u[1]
    v[2] = (u[2] - v[1]) * Int128(int_inverse_mod(m[1], m[2])) % m[2]
    return v[1] + v[2]*m[1]
end 