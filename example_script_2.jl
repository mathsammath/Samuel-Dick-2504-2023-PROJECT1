# Illustrate key functionality of the operations available via the software.

using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

# Constructing Polynomials 
println("Suppose we want to construct a polynomial. To do so, we can first construct
the polynomial x.")
@show x = x_poly()
println("We can construct larger polynomials using this x polynomial.")
@show p1 = x^2 + x + 1
@show p2 = x^7 + 2

println()

# Basic Polynomial Operations 
println("We can perform basic polynomial operations on p1 and p2 such as: addition,
subtraction, multiplication & division.")
@show +(p1, p2)
@show -(p1, p2)
@show *(p1, p2)
@show ÷(p1, p2)

println()
    
# More advanced Polynomial operations 
println("In addition to this, we can also perform more advanced operations such as:
taking the gcd of polynomials, modular arithmetic over polynomials & computing derivatives.")
@show gcd(p1, p2, 3)
@show mod(p1, 7)
@show pow_mod(p1, 4, 13)
@show derivative(p1)

println()

println("Factorization of polynomials is a key part of the code. We can 
factor polynomials mod primes & expand factorizations")
@show factorization = factor(p1, 3)
@show expand_factorization(factorization)

println()

println("Finally, we have the implementation of some basic algorithms on 
integers. Consider then, the integers 2, 5 and 6.")

println("We can take the quotient of 2 and 6:")
@show quo(2,6)

println("We can apply the Euclidean algorithm on all three numbers:")
@show euclid_alg(2,5,6)

println("We could also take the integer inverse symmetric mod of 2 and 6:")
@show int_inverse_mod(2, 6)





