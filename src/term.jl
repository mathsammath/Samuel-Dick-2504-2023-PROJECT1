#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

using Subscripts #Add this package to do superscripts.

"""
A term.
"""
struct Term  #structs are immutable by default
    coeff::Int
    degree::Int
    function Term(coeff::Int, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

"""
Creates the zero term.
"""
zero(::Type{Term})::Term = Term(0,0)

"""
Creates the unit term.
"""
one(::Type{Term})::Term = Term(1,0)

###########
# Display #
###########

"""
Show a term.
"""
#Updated this function to improve "pretty printing" of terms.
#e.g. 2x^0 → 2, 3x^1 → 3x, 1x^4 → x⁴ 
function show(io::IO, t::Term)
    if t.degree == 0 #constant 
        print(io, t.coeff)
    elseif t.degree == 1 #x ? -x : ax¹
        t.coeff == 1 ? print(io, "x") : t.coeff == - 1 ? print(io, "-x") : print(io, "$(t.coeff)⋅x", super("$(t.degree)"))
    elseif abs(t.coeff) == 1#xᵃ or -x^å
        t.coeff == 1 ? print(io, "x", super("$(t.degree)")) : print(io, "-x", super("$(t.degree)"))#-xᵃ, a ≠ 1.
    else
        print(io, "$(t.coeff)⋅x", super("$(t.degree)")) #all other cases.
    end 
end 

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms.
"""
isless(t1::Term,t2::Term)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a term at a point x.
"""
evaluate(t::Term, x::T) where T <: Number = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term,t2::Term)::Term
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term,) = Term(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::Term, t2::Term)::Term = t1 + (-t2) 

"""
Multiply two terms.
"""
*(t1::Term, t2::Term)::Term = Term(t1.coeff * t2.coeff, t1.degree + t2.degree)


"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::Term, p::Int) = Term(mod(t.coeff,p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::Term) = Term(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term,t2::Term) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term = Term(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term, n::Int) = t ÷ Term(n,0)
