#############################################################################
#############################################################################
#
# This file defines the sparse polynomial type that stores only nonzero terms 
# using a non-trivial data structure involving a dictionary and a linked list.
# Sparse polynomial type has identically functionality to the original 
# polynomial type.
#                                                                               
#############################################################################
#############################################################################

using DataStructures

####################################
# Polynomial type and construction #
####################################

"""
A Polynomial type, PolynomialSparse - designed for sparse polynomials (such as xⁿ - 1). 
"""
struct PolynomialSparse
    lst::MutableLinkedList{Term}
    dict::Dict{Int, DataStructures.ListNode{Term}}

    #Inner constructor for 0 polynomial 
    function PolynomialSparse() 
        lst = MutableLinkedList{Term}()
        append!(lst, Term(0,0)) #Append 0 term.
        dict = Dict{Int, DataStructures.ListNode{Term}}(0 => lst.node.next)
        return new(lst, dict)
    end

    #Inner constructor that creates a sparse polynomial based on arbitrary list of terms
    function PolynomialSparse(vt::Vector{Term})
        lst = MutableLinkedList{Term}()
        dict = Dict{Int, DataStructures.ListNode{Term}}()
        if isempty(vt)
            vt = [zero(Term)]
        end
        for t in vt
            insert_sorted!(lst, dict, t.degree, t)
        end 
        return new(lst, dict)
    end     
end 

#=
"""
This function maintains the invariant of the PolynomialSparse type so that 
there are no zero terms beyond the highest non-zero term.
"""
function trim!(p::PolynomialSparse)::PolynomialSparse
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end
=#

"""
Construct a sparse polynomial with a single term.
"""
PolynomialSparse(t::Term) = PolynomialSparse([t])

"""
Construct a sparse polynomial of the form x^p-x.
"""
cyclotonic_polynomial_sparse(p::Int) = PolynomialSparse([Term(1,p), Term(-1,0)])

"""
Construct a sparse polynomial of the form x-n.
"""
linear_monic_polynomial_sparse(n::Int) = PolynomialSparse([Term(1,1), Term(-n,0)])

"""
Construct a sparse polynomial of the form x.
"""
x_poly_sparse() = PolynomialSparse(Term(1,1))

"""
Creates the zero sparse polynomial.
"""
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()

"""
Creates the unit sparse polynomial.
"""
one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
one(p::PolynomialSparse) = one(typeof(p))

"""
Generates a random sparse polynomial.
"""
function rand(::Type{PolynomialSparse} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialSparse( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a sparse polynomial. May need to change this to better suit sparse. 
"""
#Updated to improve "pretty printing" of polynomials.
#e.g. 3 + 4x^2 + -3x^3 → -3x³ + 4x² + 3
function show(io::IO, p::PolynomialSparse) 
    if iszero(p)
        print(io,"0")
    else
        terms = sort([get_element(p.lst, p.dict, i) for i in keys(p.dict)]) #repeatability of code. create array of terms.
        count = 0 #iterate through terms for now...
        lowest_to_highest == false ? ordering = reverse(terms) : ordering = terms #ordering of terms depends on state of lowest_to_highest variable.
        for (i,t) in enumerate(ordering)
            if !iszero(t)
                if count == 0 #first term
                    print(io, t)
                else #all other terms
                    sign(t.coeff) > 0 && print(io, " + ", Term(abs(t.coeff), t.degree)) #Terms with positive coefficients are "added"
                    sign(t.coeff) < 0 && print(io, " - ", Term(abs(t.coeff), t.degree)) #Terms with negative coefficients are "subtracted"
                end 
                count += 1
            end
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::Polynomial, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse) = length(p.lst) 

"""
The leading term of the polynomial.
"""
leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : last(p.terms) 

"""
Returns the coefficients of the polynomial. (lowest to highest)
"""
coeffs(p::PolynomialSparse)::Vector{Int} =sort([get_element(p.lst, p.dict, i).coeff for i in keys(p.dict)])

"""
The degree of the polynomial.
"""
degree(p::Polynomial)::Int = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Polynomial)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Polynomial, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

##################################################

"""
Check if the polynomial is zero. FIX THIS
"""
iszero(p::PolynomialSparse)::Bool = p.dict[0] == [Term(0,0)]


##########################################################
##########################################################
#
# Helper functions for building the PolynomialSparse type
#                                                                               
##########################################################
##########################################################

##################################################
# Enforcing sorting and connecting to dictionary #
##################################################

"""
Assumes that the element `V` has an order (i.e. you can use <).

Assumes that the linked list `lst` is already sorted.

Assumes the dictionary `dict` does not have the key `key`.

Inserts the new element `value` to the list in its sorted position and updates the dictionary with `key` to point at the node of the linked list.
"""
function insert_sorted!(    lst::MutableLinkedList{V}, 
                            dict::Dict{K, DataStructures.ListNode{V}},
                            key::K,
                            value::V)::Nothing where {K,V}

    #Note that MutableLinkedList is implemented as a doubly pointed linked list 
    #The element lst.node is a root which points at the head and the tail.
    # lst.node.prev is the last element of the list
    # lst.node.next is the first element of the list

    haskey(dict, key) && error("Key is already in dict")
    
    #If list is empty or the value is greater than end of list, push at end
    if isempty(lst) || last(lst) <= value
        push!(lst, value)
        dict[key] = lst.node.prev #point to last since value just added to last
        return nothing
    end
    
    #if here then lst is not empty
    current_node = lst.node.next #point at first node
    i = 1
    #iterate to find 
    while current_node.data <= value
        # if current_node == lst.node.prev #if got to last node
        #     push!(lst, value) #just push at end
        #     dict[key] = lst.node.prev #point to last
        #     return nothing
        # end
        current_node = current_node.next #move to next node
    end

    #if here then current_node points at right place
    new_node = DataStructures.ListNode{V}(value) #create a new node

    #tie new_node between current_node.prev and current_node
    new_node.prev = current_node.prev 
    new_node.next = current_node

    #tie prev to new_node
    new_node.prev.next = new_node

    #tie next to new_node
    new_node.next.prev = new_node

    lst.len += 1
    dict[key] = new_node
    return nothing
 end

"""
Returns the value associated with the key `key` or `nothing` if not in list.
"""
get_element(    lst::MutableLinkedList{V}, 
                dict::Dict{K, DataStructures.ListNode{V}},
                key::K) where {K,V} = haskey(dict, key) ? dict[key].data : nothing