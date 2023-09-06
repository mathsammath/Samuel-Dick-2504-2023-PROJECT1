# Illustrate key functionality of the operations available via the software.

using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

println("The code implements same basic algorithms on integers. For example, consider the integers, 2, 5 and 6. \n")

println("We can take the quotient of 2 and 6:")
@show quo(2,6)

println()

println("We can apply the Euclidean algorithm on all three numbers:")
@show euclid_alg(2,5,6)

println()

println("We could also take the integer inverse symmetric mod of 2 and 6:")
@show int_inverse_mod(2, 6)

println()

println("We can construct polynomials. Consider the following prime polynomials denoted by p1 and p2. \n")
x = x_poly()
@show p1 = x^3 + x + 1
@show p2 = x^23 + 2

println("We can add polynomials, perform polynomial division, take the greatest-common-divisor (gcd) 
and multiply polynomials together. \n")

@show +(p1, p2)

@show ÷(p1, p2)

@show gcd(p1, p2, 3)

@show *(p1, p2)

println() 

println("In addition to the functionality detailed above, we can also perform scalar multiplication, 
check if two polynomials are the same and take the mod of a polynomial with an integer. \n")

@show 2*p1

@show ==(p1, p2)

@show mod(p1, 3)

println() 

println("Another important piece of functionality is the ability to factor polynomials. \n")
@show factor(p1, 3)


