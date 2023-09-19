####################################################################
####################################################################
#
# This file implements polynomial division for PolynomialModP types
#                                                                               
####################################################################
####################################################################

"""  
Modular algorithm for PolynomialModP types.
Assumes polynomials are over the same field. 
"""
function divide(num::PolynomialModP, den::PolynomialModP)
    return divide(num.s_poly, den.s_poly) 
end

"""
The quotient from polynomial division of two PolynomialModP types.
Assumes polynomials are over the same field. 
"""
function ÷(num::PolynomialModP, den::PolynomialModP)::PolynomialModP 
    return PolynomialModP(first(divide(num, den)(num.prime_mod)), num.prime_mod)
end 

"""
The remainder from (sparse) polynomial division. Returns a function of an integer.
Assumes polynomials are over the same field. 
"""
function rem(num::PolynomialModP, den::PolynomialModP)::PolynomialModP
    return PolynomialModP(last(divide(num,den)(num.prime_mod)), num.prime_mod)
end 