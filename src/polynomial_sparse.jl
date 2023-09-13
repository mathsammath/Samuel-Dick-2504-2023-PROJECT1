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
        #append!(lst, Term(0,0)) #Append 0 term.
        dict = Dict{Int, DataStructures.ListNode{Term}}() #(0 => lst.node.next)
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

"""
This function maintains the invariant of the Sparse Polynomial type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::PolynomialSparse)::PolynomialSparse
    i = length(p.lst)
    while i > 1
        if iszero(p.lst[i])
            pop!(p.lst)
        else
            break
        end
        i -= 1
    end
    return p
end


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
#unsure what this is actually doing???!!!
iterate(p::PolynomialSparse, state=1) = iterate([i for i in p.lst], state)

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
leading(p::PolynomialSparse)::Term = isempty(p.lst) ? zero(Term) : last(p.lst) 

"""
Returns the coefficients of the polynomial. (lowest to highest)
"""
coeffs(p::PolynomialSparse)::Vector{Int} = [i.coeff for i in p.lst]

"""
The degree of the polynomial.
"""
degree(p::PolynomialSparse)::Int = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparse)::Int = euclid_alg([i.coeff for i in p.lst])

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialSparse, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the sparse polynomial.
"""
#If term of same degree as new term exists in polynomial, throw error otherwise, add term to polynomial. 
function push!(p::PolynomialSparse, t::Term) 
    get_element(p.lst, p.dict, t.degree) === nothing ? insert_sorted!(p.lst, p.dict, t.degree, t) : 
        throw(ErrorException("Term with degree $(t.degree) already in polynomial."))
end 

"""
Pop the leading term out of the sparse polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialSparse)::Term 
    popped_term = leading(p) #popped term in leading term of polynomial
    if iszero(p) #if polynomial is zero polynomial, it must remain zero polynomial
        push!(p, zero(Term))
    end 
    delete_element!(p.lst, p.dict, leading(p).degree) #helper function to delete element 
    return popped_term
end

"""
Check if the polynomial is zero. 
"""
iszero(p::PolynomialSparse)::Bool = p.lst == zero(PolynomialSparse).lst

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a sparse polynomial.
"""
-(p::PolynomialSparse) = Polynomial(map((pt)->-pt, p.lst))

"""
Create a new sparse polynomial which is the derivative of the sparse polynomial.
"""
function derivative(p::PolynomialSparse)::PolynomialSparse 
    der_p = PolynomialSparse()
    for term in p.lst
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

#= Need to do operators before doing this!
"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialSparse) = p ÷ content(p)

"""
A square free polynomial.
"""
square_free(p::Polynomial, prime::Int)::Polynomial = (p ÷ gcd(p,derivative(p),prime))(prime)
=#

#################################
# Queries about two polynomials #
#################################

"""
Check if two sparse polynomials are the same
"""
==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.lst == p2.lst

"""
Check if a sparse polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialSparse, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two sparse polynomials.
"""
-(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse = p1 + (-p2)

"""
Multiplication of sparse polynomial and term.
"""
*(t::Term, p1::PolynomialSparse)::PolynomialSparse = iszero(t) ? PolynomialSparse() : PolynomialSparse([i for i in map((pt)->t*pt, p1.lst)])
*(p1::PolynomialSparse, t::Term)::PolynomialSparse = t*p1

"""
Multiplication of a sparse polynomial and an integer.
"""
*(n::Int, p::PolynomialSparse)::PolynomialSparse = p*Term(n,0)
*(p::PolynomialSparse, n::Int)::PolynomialSparse = n*p

"""
Integer division of a sparse polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialSparse, n::Int) = (prime)->Polynomial(map((pt)->((pt ÷ n)(prime)), p.lst)) #not sure if this is even working???????
"""
Take the mod of a sparse polynomial with an integer.
"""
function mod(f::PolynomialSparse, p::Int)::PolynomialSparse
    f_out = PolynomialSparse()
    for i in f.lst
        f_out += mod(i,p) 
    end 
    return f_out        
end

"""
Power of a sparse polynomial mod prime.
"""
function pow_mod(p::PolynomialSparse, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end

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

"""
Assumes the dictionary `dict` has the key `key`. and it is pointing at the linked list.
"""
function delete_element!(   lst::MutableLinkedList{V}, 
                            dict::Dict{K, DataStructures.ListNode{V}},
                            key::K)::Nothing where {K,V}
    haskey(dict, key) || error("Key is not in dict")

    node = dict[key]
    delete!(dict, key)

    node.prev.next = node.next
    node.next.prev = node.prev
    lst.len -= 1

    return nothing
end
