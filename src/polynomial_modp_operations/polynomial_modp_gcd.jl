#############################################################################
#############################################################################
#
# This file implements polynomial GCD for polynomials of type PolynomialModP
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for PolynomialModP types over p.
Assumes polynomials are over the same field.
"""
#not 100% sure if this is working!!!!
function extended_euclid_alg(a::PolynomialModP, b::PolynomialModP) #remove prime
    s_type = typeof(a.s_poly) #type either sparse of sparse128
    old_r, r = a, b
    old_s, s = one_modP(s_type, a.prime_mod), zero_modP(s_type, a.prime_mod)
    old_t, t = zero_modP(s_type, a.prime_mod), one_modP(s_type, a.prime_mod)

    while !iszero(r)
        q = ÷(old_r, r)
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    end
    g, s, t = old_r, old_s, old_t

    @assert s*a + t*b - g == 0
    return g, s, t  
end

"""
The GCD of two polynomials of type PolynomialModP modulo p.
Assumes polynomials are over the same field.
"""
gcd(a::PolynomialSparse, b::PolynomialSparse) = extended_euclid_alg(a,b) |> first