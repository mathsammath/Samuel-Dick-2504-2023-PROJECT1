#####################################################################################
#####################################################################################
#
# This file contains units tests for polynomial operations for sparse polynomial type
# Randomly generated polynomials can be dense, i.e. x⁷ + 6x⁶ + 4x⁵ + 3x⁴ + 2x²
# As such, the number of iterations of some tests relative to dense 
# testing has been reduced.
#                                                                               
#####################################################################################
#####################################################################################


"""
Test product of sparse polynomials.
"""
function prod_test_poly_sparse(;N::Int = 30, N_prods::Int = 20, seed::Int = 0) #obviously taking a very long time (relatively)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialSparse(Term(1,0))
        for _ in 1:N_prods
            p = rand(PolynomialSparse)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly - PASSED")
end

"""
Test derivative of sparse polynomials (as well as product).
"""
function prod_derivative_test_poly_sparse(;N::Int = 30,  seed::Int = 0)
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
function division_test_poly_sparse(;prime::Int = 101, N::Int = 10^3, seed::Int = 0) 
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
function ext_euclid_test_poly_sparse(;prime::Int=101, N::Int = 10, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly - PASSED")
end