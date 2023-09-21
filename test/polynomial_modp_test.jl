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
function prod_test_poly_modp(;N::Int = 10^4, N_prods::Int = 20, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)

        prime_num = rand(primes(1,20)) #random prime p in range 1:20
        f1 = PolynomialModP(p1, prime_num) #initialise modp polynomial over field Zₚ
        f2 = PolynomialModP(p2, prime_num) #as above

        prod_sparse = p1*p2
        prod = f1*f2 
        #Product of PolynomialModP polynomials should be equivalent to product of 
        #sparse polynomials, mod the field which polynomials are working over.
        @assert prod.s_poly.lst == mod(prod_sparse, prime_num).lst
    end
    println("prod_test_poly_modp - PASSED")
end

"""
Test derivative of PolynomialModP polynomials (as well as product).
"""
function prod_derivative_test_poly_modp(;N::Int = 10^3,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly_modp - PASSED")
end

"""
Test ÷ and "rem" functions for PolynomialModP types.
"""
function div_rem_test_poly_modp(;prime::Int = 101, N::Int = 10^4, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)

        prime_num = rand(primes(1,20)) #random prime p in range 1:20
        f1 = PolynomialModP(p1, prime_num) #initialise modp polynomial over field Zₚ
        f2 = PolynomialModP(p2, prime_num) #as above

        f_prod = f1*f2

        q, r = PolynomialModP(PolynomialSparse, prime_num), PolynomialModP(PolynomialSparse, prime_num)
        try
            #cases where one of q or r is a constant (consequence of rand function)
            leading(q.s_poly).degree == 0 || leading(r.s_poly).degree == 0 ? continue : nothing 
            q, r = ÷(f_prod,f2), rem(f_prod,f2) #equivalent to result of div(q,r)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert f2 == 0
            else
                throw(e)
            end
        end
        @assert iszero(q*f2+r - f_prod)
    end
    println("division_test_poly_modp - PASSED")
end

"""
Test the push! method for PolynomialModP
"""
function push_test_poly_modp(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialModP)
        coeff = rand(-50:50)
        deg = rand(0:25)
        t = Term(coeff, deg)
        try
            push!(p, t) #add the term to the polynomial
            @assert t ∈ p.s_poly.lst
        catch e #if term of same degree already exists, throw error. 
            if typeof(e) == ErrorException
                @assert true 
            end 
        end
    end 
    println("push_test_poly_modp - PASSED")
end 


"""
Test the pop! method for sparse polynomials 
"""
function pop_test_poly_modp(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialModP)
        lt = leading(p)
        pop!(p)
        @assert lt ∉ p.s_poly.lst #assert leading term has been removed from p1
    end 
    println("pop_test_poly_modp - PASSED")
end 