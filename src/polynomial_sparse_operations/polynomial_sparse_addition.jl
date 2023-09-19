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
    iszero(t) && return p
    p = deepcopy(p)
    if t.degree ∈ keys(p.dict) #polynomial contains term of same degree
        add_term = get_element(p.lst, p.dict, t.degree) + t #sum terms of same degree 
        delete_element!(p.lst, p.dict, t.degree) 
        push!(p, add_term) #add new term in place of removed term 
    else 
        push!(p, t)
    end 
    return p
end

+(t::Term, p::PolynomialSparse) = p + t

"""
Add two polynomials of type PolynomialSparse.
"""
function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    q1 = deepcopy(p1) #pop! and push! mutate polynomials
    q2 = deepcopy(p2)
    sum_pol = PolynomialSparse() #initialise sum of poylnomials 
    while !iszero(q1) || !iszero(q2)
        if degree(q1) == degree(q2) #degrees equal, add leading terms
            ltq1, ltq2 = pop!(q1), pop!(q2)
            if !iszero(ltq1 + ltq2) #PolynomialSparse contains only nonzero terms
                push!(sum_pol, ltq1 + ltq2)
            end 
        elseif degree(q1) > degree(q2)  
            push!(sum_pol, pop!(q1)) #leading term of q1 added  
        else 
            push!(sum_pol, pop!(q2)) #leading term of q2 added
        end 
    end 
    return sum_pol 
end

"""
Add a sparse polynomial and an integer.
"""
+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)