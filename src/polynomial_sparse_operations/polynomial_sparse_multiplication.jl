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
    p_out = PolynomialSparse()
    for t in p1.lst
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a sparse polynomial.
"""
function ^(p::PolynomialSparse, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end
