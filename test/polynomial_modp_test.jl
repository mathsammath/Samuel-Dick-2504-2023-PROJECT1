#####################################################################################
#####################################################################################
#
# This file contains units tests for polynomial operations for PolynomialModP type
#                                                                               
#####################################################################################
#####################################################################################
using Primes 

"""
Test product of polynomials of type PolynomialModP.
"""
function prod_test_poly_modp(;N::Int = 30, N_prods::Int = 20, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        println(p1)
        println(p2)

        prime_num = rand(primes(1,20)) #random prime p in range 1:20
        f1 = PolynomialModP(p1, prime_num) #initialise modp polynomial over field Zₚ
        f2 = PolynomialModP(p2, prime_num) #initialise modp polynomial over field Zₚ
        println(f1)
        println(f2)

        prod_sparse = p1*p2
        prod = f1*f2 
        #Product of PolynomialModP polynomials should be equivalent to product of 
        #sparse polynomials, mod the field which polynomials are working over.
        prod == mod(prod_sparse, prime_num)
    end
    println("prod_test_poly_modp - PASSED")
end

#=
"""
Test derivative of sparse polynomials (as well as product).
"""
function prod_derivative_test_poly_modp(;N::Int = 30,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly - PASSED")
end


"""
Test division of sparse polynomials modulo p.
"""
function division_test_poly_modp(;prime::Int = 101, N::Int = 10^3, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p_prod = p1*p2
        q, r = PolynomialSparse(), PolynomialSparse()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("division_test_poly - PASSED")
end

"""
Test the extended euclid algorithm for sparse polynomials modulo p.
"""
function ext_euclid_test_poly_modp(;prime::Int=101, N::Int = 10, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly - PASSED")
end

"""
Test the push! method for sparse polynomials
"""
function push_test_poly_modp(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparse)
        coeff = rand(-50:50)
        deg = rand(0:25)
        t = Term(coeff, deg)
        try
            push!(p, t) #add the term to the polynomial
            @assert t ∈ p.lst
        catch e #if term of same degree already exists, throw error. 
            if typeof(e) == ErrorException
                @assert true 
            end 
        end
    end 
    println("push_test_poly_sparse - PASSED")
end 


"""
Test the pop! method for sparse polynomials 
"""
function pop_test_poly_modp(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparse)
        lt = leading(p)
        pop!(p)
        @assert lt ∉ p.lst #assert leading term has been removed from p1
    end 
    println("pop_test_poly_sparse - PASSED")
end 
=#