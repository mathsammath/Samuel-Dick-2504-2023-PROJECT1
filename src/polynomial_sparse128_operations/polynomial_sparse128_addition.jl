#######################################################################
#######################################################################
#
# This file implements polynomial addition for PolynomialSparse128 types.
#                                                                               
########################################################################
########################################################################

"""
Add a sparse polynomial and a term.
"""
function +(p::PolynomialSparse128, t::Term128)
    p = deepcopy(p)
    if t.degree ∈ [i.degree for i in p.lst] #Term of same degree exists in p
         for i in 1:length(p.lst)
            if p.lst[i].degree == t.degree
                p.lst[i] += t
            end 
        end 
    else #Term of same degree does not exist in p
        push!(p, t)
    end
    return trim!(p)
end

+(t::Term, p::PolynomialSparse128) = p + t

"""
Add two sparse polynomials.
"""
function +(p1::PolynomialSparse128, p2::PolynomialSparse128)::PolynomialSparse128
    p = deepcopy(p1)
    for t in p2.lst 
        p += t
    end
    return p
end

"""
Add a sparse polynomial and an integer.
"""
+(p::PolynomialSparse128, n::Int) = p + Term128(Int128(n),0)
+(n::Int, p::PolynomialSparse128) = p + Term128(Int128(n),0)