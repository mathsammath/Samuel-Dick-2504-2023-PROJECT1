#############################################################################
#############################################################################
#
# This file implements polynomial GCD for dense polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials modulo prime.
"""
function extended_euclid_alg(a::PolynomialDense, b::PolynomialDense, prime::Int)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialDense), zero(PolynomialDense)
    old_t, t = zero(PolynomialDense), one(PolynomialDense)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::PolynomialDense, b::PolynomialDense, prime::Int) = extended_euclid_alg(a,b,prime) |> first