#######################################################################
#######################################################################
#
# This file implements polynomial addition for PoylnomialModP types.
#                                                                               
########################################################################
########################################################################

"""
Add a PolynomialModP and a term.
"""
function +(p::PolynomialModP, t::Term)
    return PolynomialModP(p.s_poly + t , p.prime_mod)
end

+(t::Term, p::PolynomialSparse) = p + t

"""
Add two polynomials of type PolynomialModP.
Assumes polynomials are over the same field.
"""
function +(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    return PolynomialModP(p1.s_poly + p2.s_poly, p1.prime_mod)
end

"""
Add a PolynomialModP and an integer.
"""
+(p::PolynomialModP, n::Int) = PolynomialModP(p.s_poly + n, p.prime_mod)
+(n::Int, p::PolynomialModP) = PolynomialModP(p.s_poly + n, p.prime_mod)