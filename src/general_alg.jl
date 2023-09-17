#############################################################################
#############################################################################
#
# This file implement several basic algorithms for integers. 
#                                                                               
#############################################################################
#############################################################################

"""
The quotient between two numbers
"""
function quo(a::T,b::T) where T <: Integer
    a < 0 && return -quo(-a,b)
    a < b && return 0
    return 1 + quo(a-b, b)
end

"""
Euclid's algorithm.
"""
function euclid_alg(a, b; rem_function = %)
    b == 0 && return a
    return euclid_alg(b, rem_function(a,b); rem_function = rem_function)
end

"""
Euclid's algorithm on multiple arguments.
"""
euclid_alg(a...) = foldl((a,b)->euclid_alg(a,b), a; init = 0)

"""
Euclid's algorithm on a vector.
"""
euclid_alg(a::Vector{T}) where T <: Integer = euclid_alg(a...)


"""
The extended Euclidean algorithm.
"""
function ext_euclid_alg(a, b, rem_function = %, div_function = ÷)
    a == 0 && return b, 0, 1
    g, t, s = ext_euclid_alg(rem_function(b,a), a, rem_function, div_function)
    s = s - div_function(b,a)*t
    @assert g == a*s + b*t
    return g, s, t
end

"""
Display the result of the extended Euclidean algorithm.
"""
pretty_print_egcd((a,b),(g,s,t)) = println("$a × $s + $b × $t = $g") #\times + [TAB]

"""
Integer inverse symmetric mod
"""
function int_inverse_mod(a::Int, m::Int)::Int 
    if mod(a, m) == 0
        error("Can't find inverse of $a mod $m because $m divides $a") 
    end
    return mod(ext_euclid_alg(a,m)[2],m)
end