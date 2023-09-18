#############################################################################
#############################################################################
#
# This file defines the PolynomialModP type
#                                                                               
#############################################################################
#############################################################################

using DataStructures, Primes

####################################
# Polynomial type and construction #
####################################

"""
A Polynomial type, PolynomialModP- designed for sparse polynomials with the 
addition of working with coefficient arithmetic over Zₚ.  
"""
struct PolynomialModP
    prime_mod::Int #prime
    s_poly::Union{PolynomialSparse, PolynomialSparse128} #sparse polynomial 

    #Inner constructor for 0 PolynomialSparse
    function PolynomialModP(::Type{PolynomialSparse}, p::Int)
        prime_mod = p
        s_poly = PolynomialSparse()
        return new(prime_mod, s_poly)
    end 

    #Inner constructor for 0 PolynomialSparse128
    function PolynomialModP(::Type{PolynomialSparse128}, p::Int)
        prime_mod = p
        s_poly = PolynomialSparse128()
        return new(prime_mod, s_poly)
     end 

    #Inner constructor for PolynomialModP based on arbitrary list of terms (Int) and a prime
    function PolynomialModP(vt::Vector{Term}, p::Int)
        prime_mod = p
        s_poly = mod(PolynomialSparse(vt), prime_mod) #ensure polynomial is constructed with coefficients mod p
        return new(prime_mod, s_poly)
    end  

    #Inner constructor for PolynomialModP based on arbitrary list of terms (Int128) and a prime
    function PolynomialModP(vt::Vector{Term128}, p::Int)
        prime_mod = p
        s_poly = mod(PolynomialSparse128(vt), prime_mod) #ensure polynomial is constructed with coefficients mod p
        return new(prime_mod, s_poly)
    end 
end 

"""
Construct a PolynomialModP with a single term and prime.
"""
PolynomialModP(t::Union{Term, Term128}, p::Int) = PolynomialSparse([t], p)

#=
"""
Construct a PolynomialModP of the form x^p-x given a prime.
"""
cyclotonic_polynomial_modp(p::Int, n::Int) = PolynomialSparse([Term(1,p), Term(-1,0)])

"""
Construct a PolynomialModP of the form x-n.
"""
linear_monic_polynomial_sparse(n::Int) = PolynomialSparse([Term(1,1), Term(-n,0)])
=#

"""
Construct a PolynomialModP of the form x given a prime p.
"""
x_poly_modp(::Type{PolynomialSparse}, p::Int) = PolynomialModP([Term(1,1)], p)
x_poly_modp(::Type{PolynomialSparse128}) = PolynomialModP([Term128(Int128(1),1)], p)


"""
Creates the zero PolynomialModP given some prime.
"""
zero_modP(::Type{PolynomialSparse}, p::Int)::PolynomialModP = PolynomialModP(PolynomialSparse, p)
zero_modP(::Type{PolynomialSparse128}, p::Int)::PolynomialModP = PolynomialModP(PolynomialSparse128, p)

"""
Creates the unit PolynomialModP given some prime.
"""
one_modP(::Type{PolynomialSparse}, p)::PolynomialModP = PolynomialModP(one(Term), p)
one_modP(::Type{PolynomialSparse128}, p)::PolynomialModP = PolynomialModP(one(Term128), p)

#one(p::PolynomialModP) = one(typeof(p), p.prime_mod) something funny going on.

"""
Generates a random PolynomialModP.
"""
function rand(::Type{PolynomialModP} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        prime_mod = rand(primes(1,20)) #random prime p in range [1, 100]
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialModP( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)], prime_mod )
        condition(p) && return p
    end
end

###########
# Display #
###########
"""
Show a sparse polynomial. May need to change this to better suit sparse. 
"""
function show(io::IO, p::PolynomialModP) 
    if iszero(p)
        print(io,"0")
    else
        first_term = true #first nonzero term printed differently to other terms.
        lowest_to_highest == false ? ordering = reverse(p.s_poly.lst) : ordering = p.s_poly.lst #ordering of terms depends on state of lowest_to_highest variable.
        for (i,t) in enumerate(ordering)
            if !iszero(t)
                if first_term == true 
                    print(io, t) #will print negative with sign and positive without.
                else 
                    sign(t.coeff) > 0 && print(io, " + ", Term(abs(t.coeff), t.degree)) #Terms with positive coefficients are "added"
                    sign(t.coeff) < 0 && print(io, " - ", Term(abs(t.coeff), t.degree)) #Terms with negative coefficients are "subtracted"
                end 
                first_term = false 
            end
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the PolynomialModP. This implements the iteration interface.
"""
iterate(p::PolynomialModP, state=1) = iterate(collect(p.s_poly.lst), state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the PolynomialModP.
"""
length(p::PolynomialModP) = length(p.s_poly.lst) 

"""
The leading term of the PolynomialModP.
"""
leading(p::PolynomialModP) = isempty(p.s_poly.lst) ? zero(Term) : last(p.s_poly.lst) #might need to return Term128... but maybe not??? 

"""
Returns the coefficients of the PolynomialModP. 
"""
function coeffs(p::PolynomialModP)::Vector{Int}
    coeff_v = Int[] #initialise 
    f = deepcopy(p.s_poly) #push! mutates p.
    while length(f) > 0 
        append!(coeff_v, pop!(f).coeff) 
    end 
    return coeff_v 
end 

"""
The degree of the PolynomialModP.
"""
degree(p::PolynomialModP)::Int = leading(p).degree 

"""
The content of the PolynomialModP is the GCD of its coefficients.
"""
content(p::PolynomialModP)::Int = euclid_alg(coeffs(p))

"""
Evaluate the PolynomialModP at a point `x`.
"""
function evaluate(f::PolynomialModP, x::T) where T <: Number
    evaluate(f.s_poly, x)
end 

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the PolynomialModP.
"""
#If term of same degree as new term exists in polynomial, throw error otherwise, add term to polynomial. 
function push!(p::PolynomialModP, t::Union{Term, Term128}) 
    t = mod(t, p.prime_mod) #make sure the term we're pushing is modulo p
    push!(p.s_poly, t)
end 

"""
Pop the leading term out of the PolynomialModP. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialModP)::Union{Term, Term128} 
    pop!(p.s_poly)
end

"""
Check if the PolynomialModP is zero. 
"""
iszero(p::PolynomialModP)::Bool = p.s_poly.lst == zero(PolynomialSparse).lst || p.s_poly.lst == zero(PolynomialSparse128).lst ||
                                    p.s_poly.lst == MutableLinkedList{Term}(Term(0,0)) || p.s_poly.lst == MutableLinkedList{Term128}(Term128(Int128(0),0))

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a PolynomialModP.
"""
-(p::PolynomialModP) = PolynomialModP(map(t->-t, p.s_poly), p.prime_mod)

"""
Create PolynomialModP which is the derivative of a PolynomialModP.
"""
function derivative(p::PolynomialModP)::PolynomialModP
    derivative(p.s_poly)
end

"""
The prim part (multiply a PolynomialModP by the inverse of its content).
"""
prim_part(p::PolynomialModP) = prim_part(p.s_poly)

"""
A square free PolynomialModP.
"""
#square_free(p::PolynomialSparse, prime::Int)::PolynomialSparse = (p ÷ gcd(p,derivative(p),prime))(prime) comeback to this. important for later

#################################
# Queries about two polynomials #
#################################

"""
Check if two PolynomialModP's are the same
"""
==(p1::PolynomialModP, p2::PolynomialModP)::Bool = p1.s_poly == p2.s_poly

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
function *(t::Term, p1::PolynomialSparse)::PolynomialSparse 
     iszero(t) ? PolynomialSparse() : 
     p2 = deepcopy(p1) #pop mutates p1
     p = PolynomialSparse() 
     for _ in 1:length(p2)
        push!(p, pop!(p2)*t) #multiplication done term-by-term
     end 
     return p 
end 

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
÷(p::PolynomialSparse, n::Int) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.lst)) 
"""
Take the mod of a sparse polynomial with an integer.
"""
function mod(f::PolynomialSparse, p::Int)::PolynomialSparse
    mod_p = PolynomialSparse() #iniialise 
    f_copy = deepcopy(f)
    for _ in 1:length(f)
        !iszero(mod(leading(f_copy), p)) ? push!(mod_p, mod(pop!(f_copy), p)) : pop!(f_copy)
    end 
    return mod_p
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
