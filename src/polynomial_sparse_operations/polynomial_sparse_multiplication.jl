#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for sparse polynomial types
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two sparse polynomials.
"""
function *(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p_out = PolynomialSparse() #initialise result of multiplication
    p1_terms = collect(p1.lst)
    for t in p1_terms 
        p_out += t*p2 #multiply each term of p1 by p2
    end 
    return p_out 
end

"""
Power of a sparse polynomial.
"""
function ^(p::PolynomialSparse, n::Int)
    n < 0 && error("No negative power")
    b = reverse(string(n; base=2)) #binary string representation of n, reversed.
    ans, w = Int(1), p
    for i in 1:length(b)
        if b[i] == '1'
            ans = ans*w
        end
        w = w*w
    end 
    return ans 
end
