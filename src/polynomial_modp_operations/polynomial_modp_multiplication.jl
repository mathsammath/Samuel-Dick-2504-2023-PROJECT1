#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for PolynomialModP types
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials of type PolynomialModP.
Assumes polynomials are over the same field.
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    return PolynomialModP(p1.s_poly*p2.s_poly, p1.prime_mod)
end

"""
Power of a polynomial of type PolynomialModP.
"""
function ^(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end
