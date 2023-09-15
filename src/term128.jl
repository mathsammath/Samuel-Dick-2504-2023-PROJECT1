#############################################################################
#############################################################################
#
# This file defines the Term128 type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term with integer coefficients of type Int128.
"""
struct Term128  
    coeff::Int128
    degree::Int
    function Term128(coeff::Int128, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

"""
Creates the zero term.
"""
zero(::Type{Term128})::Term128 = Term128(Int128(0),0)

"""
Creates the unit term.
"""
one(::Type{Term128})::Term128 = Term128(Int128(1),0)

###########
# Display #
###########

"""
Show a term.
"""
#Updated this function to improve "pretty printing" of terms.
#e.g. 2x^0 → 2, 3x^1 → 3x, 1x^4 → x⁴ 
function show(io::IO, t::Term128)
    if t.degree == 0 
        print(io, t.coeff) #constant 
    elseif t.degree == 1 
        t.coeff == 1 ? print(io, "x") : t.coeff == - 1 ? print(io, "-x") : print(io, "$(t.coeff)⋅x") #x, -x, ax
    elseif abs(t.coeff) == 1 
        t.coeff == 1 ? print(io, "x", super_index(t.degree)) : print(io, "-x", super_index(t.degree)) #xᵃ, -xᵃ
    else 
        print(io, "$(t.coeff)⋅x", super_index(t.degree)) #all other cases
    end
end 

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term128)::Bool = iszero(t.coeff)

"""
Compare two terms.
"""
isless(t1::Term128,t2::Term128)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a term at a point x.
"""
evaluate(t::Term128, x::T) where T <: Number = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term128,t2::Term128)::Term128
    @assert t1.degree == t2.degree
    Term128(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term128,) = Term128(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::Term128, t2::Term128)::Term128 = t1 + (-t2) 

"""
Multiply two terms.
"""
*(t1::Term128, t2::Term128)::Term128 = Term128(t1.coeff * t2.coeff, t1.degree + t2.degree)


"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::Term128, p::Int) = Term128(mod(t.coeff,p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::Term128) = Term128(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term128,t2::Term128) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term128 = Term128(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term128, n::Int) = t ÷ Term128(Int128(n),0)

##########################################################
##########################################################
#
# Helper functions for Term type
#                                                                               
##########################################################
##########################################################

"""
Returns an integer as a superscript. e.g. 1 → ¹.
"""
function super_index(n::Int)::String
    super_dict = Dict(
        '0' => "⁰",
        '1' => "¹",
        '2' => "²",
        '3' => "³",
        '4' => "⁴",
        '5' => "⁵",
        '6' => "⁶",
        '7' => "⁷",
        '8' => "⁸",
        '9' => "⁹"
    )
    n_string = string(n)
    super_string = ""
    for i in n_string 
        super_string *= super_dict[i]
    end 
    return super_string 
end 