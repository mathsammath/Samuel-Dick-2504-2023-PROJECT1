#######################################################################
#######################################################################
#
# This file implements polynomial addition for PolynomialSparse128 types.
#                                                                               
########################################################################
########################################################################

"""
Add a polynomial (sparse128) and a term.
"""
function +(p::PolynomialSparse128, t::Term128)
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

+(t::Term128, p::PolynomialSparse128) = p + t

"""
Add two sparse128 polynomials.
"""
function +(p1::PolynomialSparse128, p2::PolynomialSparse128)::PolynomialSparse128
    q1 = deepcopy(p1) #pop! and push! mutate polynomials
    q2 = deepcopy(p2)
    sum_pol = PolynomialSparse128() #initialise sum of poylnomials 
    while !iszero(q1) || !iszero(q2)
        if degree(q1) == degree(q2) #degrees equal, add leading terms
            ltq1, ltq2 = pop!(q1), pop!(q2)
            if ltq1 + ltq2 ≠ 0 #PolynomialSparse contains only nonzero terms
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
Add a sparse128 polynomial and an integer.
"""
+(p::PolynomialSparse128, n::Int128) = p + Term128(n,0)
+(n::Int128, p::PolynomialSparse128) = p + Term128(n,0)