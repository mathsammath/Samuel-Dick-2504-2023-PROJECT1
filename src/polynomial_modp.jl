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
A Polynomial type, PolynomialModP- designed for polynomials (PolynomialSparse/PolynomialSparse128) 
supporting coefficient arithmetic over Zₚ.  
"""
struct PolynomialModP
    prime_mod::Int #prime
    s_poly::Union{PolynomialSparse, PolynomialSparse128} #sparse polynomial 

    #Inner constructor for 0 PolynomialModP (PolynomialSparse) over some prime 
    function PolynomialModP(::Type{PolynomialSparse}, p::Int)
        prime_mod = p
        s_poly = PolynomialSparse()
        return new(prime_mod, s_poly)
    end 

    #Inner constructor for 0 PolynomialModP (PolynomialSparse128) over some prime 
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
Construct a PolynomialModP from a Sparse or Sparse128 type.
"""
PolynomialModP(p::Union{PolynomialSparse, PolynomialSparse128}, t::Int) = PolynomialModP(collect(p.lst), t)

"""
Construct a PolynomialModP with a single term and prime.
"""
PolynomialModP(t::Union{Term, Term128}, p::Int) = PolynomialModP([t], p)

"""
Construct a PolynomialModP of the form x^n-x given a prime.
"""
cyclotonic_polynomial_modp(::Type{PolynomialSparse}, p::Int, n::Int) = PolynomialModP([Term(1,n), Term(-1,0)], p)
cyclotonic_polynomial_modp(::Type{PolynomialSparse128}, p::Int, n::Int) = PolynomialModP([Term128(Int128(1),n), Term128(Int128(-1),0)], p)

"""
Construct a PolynomialModP of the form x-n.
"""
linear_monic_polynomial_modp(::Type{PolynomialSparse}, p::Int, n::Int) = PolynomialModP([Term(1,1), Term(-n,0)], p)
linear_monic_polynomial_modp(::Type{PolynomialSparse128}, p::Int, n::Int) = PolynomialModP([Term128(Int128(1),1), Term128(Int128(-n),0)], p)

"""
Construct a PolynomialModP of the form x given a prime p.
"""
x_poly_modp(::Type{PolynomialSparse}, p::Int) = PolynomialModP([Term(1,1)], p)
x_poly_modp(::Type{PolynomialSparse128}, p::Int) = PolynomialModP([Term128(Int128(1),1)], p)


"""
Creates the zero PolynomialModP given some prime.
"""
zero_modP(::Type{PolynomialSparse}, p::Int)::PolynomialModP = PolynomialModP(PolynomialSparse, p)
zero_modP(::Type{PolynomialSparse128}, p::Int)::PolynomialModP = PolynomialModP(PolynomialSparse128, p)

"""
Creates the unit PolynomialModP given some prime.
"""
one_modP(::Type{PolynomialSparse}, p)::PolynomialModP = PolynomialModP([one(Term)], p)
one_modP(::Type{PolynomialSparse128}, p)::PolynomialModP = PolynomialModP([one(Term128)], p)

one(p::PolynomialModP) = one_modP(typeof(p.s_poly), p.prime_mod) 

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
        prime_mod = rand(primes(1,20)) #random prime p in range 1:20
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
    show(io, p.s_poly)
    print(" (mod ", p.prime_mod, ")")
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
function coeffs(p::PolynomialModP)::Vector{Union{Int, Int128}}
    coeffs(p.s_poly)
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
    PolynomialModP(derivative(p.s_poly), p.prime_mod)
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
==(p::PolynomialModP, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two PolynomialModP's.
"""
-(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP = p1 + (-p2)

"""
Multiplication of a PolynomialModP and a term.
"""
function *(t::Term, p1::PolynomialModP)::PolynomialModP
    PolynomialModP(collect((*(t, p1.s_poly)).lst) , p1.prime_mod)
end 

*(p1::PolynomialModP, t::Term)::PolynomialModP = t*p1

"""
Multiplication of a PolynomialModP and an integer.
"""
function *(n::Int, p::PolynomialModP)::PolynomialModP
     typeof(p.s_poly) == PolynomialSparse ? p*Term(n,0) : p*Term128(Int128(n), 0)
end 

*(p::PolynomialModP, n::Int)::PolynomialModP = n*p

"""
Integer division of a PolynomialModP by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
#not working 
÷(p::PolynomialModP, n::Int) = ÷(p.s_poly, n)(p.prime_mod)
"""
Take the mod of a PolynomialModP with an integer.
"""
function mod(f::PolynomialModP, p::Int)::PolynomialModP
    return PolynomialModP(f.s_poly, p)
end

"""
Power of a PolynomialModP mod prime.
"""
function pow_mod(p::PolynomialModP, n::Int, prime::Int)
    x = PolynomialModP(p.s_poly, prime) #since we now want to work over a new prime 
    return ^(x, n)
end
