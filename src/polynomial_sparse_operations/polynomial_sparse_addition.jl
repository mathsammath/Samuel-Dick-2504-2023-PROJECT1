#######################################################################
#######################################################################
#
# This file implements polynomial addition for sparse polynomial types.
#                                                                               
########################################################################
########################################################################

"""
Add a sparse polynomial and a term.
"""
function +(p::PolynomialSparse, t::Term)
    p = deepcopy(p)
    if t.degree ∈ [i.degree for i in p.lst] #Term of same degree exists in p
         for i in 1:length(p.lst) #this is not correct. must change.
            if p.lst[i].degree == t.degree
                p.lst[i] += t
            end 
        end 
    else #Term of same degree does not exist in p
        push!(p, t)
    end
    return trim!(p)
end

+(t::Term, p::PolynomialSparse) = p + t

"""
Add two sparse polynomials.
"""
function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2.lst 
        p += t
    end
    return p
end

"""
Add a sparse polynomial and an integer.
"""
+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)