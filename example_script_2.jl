# Illustrate key functionality of the operations available via the software.

using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

# Constructing Polynomials 
println("Suppose we want to construct a polynomial of sparse type. To do so, we can first construct
the polynomial x.")
@show x = x_poly_sparse()
println("We can construct larger polynomials using this x polynomial. Note these polynomials 
have purposefully been defined over primes.")
@show p1 = mod(rand(PolynomialSparse), 7)
@show p2 = mod(rand(PolynomialSparse), 3)

println()

# Basic Polynomial Operations 
println("We can perform basic polynomial operations on p1 and p2 such as: addition,
subtraction, multiplication & division.")
@show +(p1, p2)
@show -(p1, p2)
@show *(p1, p2)
@show ÷(p1, p2)
println("Note that the ÷ function returns a tuple, (quotient, remainder) of the polynomial division.")

println()
     
println("In addition to this, we can also perform more advanced operations such as:
taking powers of polynomials, modular arithmetic over polynomials & computing derivatives.")
@show mod(p1, 7)
@show pow_mod(p1, 4, 13)
@show derivative(p2)

println()

println("We also have the implementation of some basic algorithms on 
integers. Consider then, the integers 2, 5 and 6.")

println("We can take the quotient of 2 and 6:")
@show quo(2,6)

println("We can apply the Euclidean algorithm on all three numbers:")
@show euclid_alg(2,5,6)

println("We could also take the integer inverse symmetric mod of 2 and 6:")
@show int_inverse_mod(2, 6)
println()

println("The file also has operations for polynomial factorisation, these will be implemented later for the PolynomialSparse type.")




