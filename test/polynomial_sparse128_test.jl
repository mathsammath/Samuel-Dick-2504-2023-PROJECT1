#####################################################################################
#####################################################################################
#
# This file contains units tests for polynomial operations for PolynomialSparse128.
# Similar testing to PolynomialSparse but edge cases investigated (relative to 
# PolynomialSparse).
#                                                                               
#####################################################################################
#####################################################################################

"""
Test product for polynomials of type PolynomialSparse128.
"""
function prod_test_poly_sparse128(;N::Int = 30, N_prods::Int = 20, seed::Int = 0) #obviously taking a very long time (relatively)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialSparse128(Term128(Int128(1),0))
        for _ in 1:N_prods
            p = rand(PolynomialSparse128)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly_sparse128 - PASSED")
end

"""
Test derivative for polynomials of type PolynomialSparse128.
"""
function prod_derivative_test_poly_sparse128(;N::Int = 30,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly_sparse128 - PASSED")
end


"""
Test division for polynomials of type PolynomialSparse128.
"""
function division_test_poly_sparse128(;prime::Int = 101, N::Int = 10^3, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        p_prod = p1*p2
        q, r = PolynomialSparse128(), PolynomialSparse128()
        try
            #cases where one of p_prod or p2 is a constant (consequence of rand function)
            leading(p_prod).degree == 0 || leading(p2).degree == 0 ? continue : nothing             
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
    println("division_test_poly_sparse128 - PASSED")
end

"""
Test the extended euclid algorithm for polynomials of type PolynomialSparse128.
"""
function ext_euclid_test_poly_sparse128(;prime::Int=101, N::Int = 30, seed::Int = 0) 
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly_sparse128 - PASSED")
end

"""
Test the push! method for polynomials of type PolynomialSparse128.
"""
function push_test_poly_sparse128(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparse128)
        coeff = Int128(rand(-50:50))
        deg = rand(0:25)
        t = Term128(coeff, deg)
        try
            push!(p, t) #add the term to the polynomial
            @assert t ∈ p.lst
        catch e #if term of same degree already exists, throw error. 
            if typeof(e) == ErrorException
                @assert true 
            end 
        end
    end 
    println("push_test_poly_sparse128 - PASSED")
end 


"""
Test the pop! method for polynomials of type PolynomialSparse128.
"""
function pop_test_poly_sparse128(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparse128)
        lt = leading(p)
        pop!(p)
        @assert lt ∉ p.lst #assert leading term has been removed from p1
    end 
    println("pop_test_poly_sparse128 - PASSED")
end 


"""
Test the mod method for polynomials of type PolynomialSparse128.
Polynomials tested mod 7 and mod 5. 
Examples hardcoded. 
"""
function mod_test_poly_sparse128(p::Int = 7, q::Int = 5)
    x = x_poly_sparse128() 
    p1 = Int128(3)x^7 + Int128(7)x^6 + Int128(10)x^5 + x^4 + Int128(4)x^3 + Int128(2)x^2
    p2 = x^7 + Int128(2)x^6 + Int128(20)x^5 + Int128(6)x^4 + Int128(2)x^3 + Int128(13)x^2
    p3 = x^7 + x^6 + Int128(7)x^5 + Int128(30)x^4 + Int128(70)x^3 + Int128(140)x^2

    #mod(p1, 7) = 3⋅x⁷ + 3⋅x⁵ + x⁴ + 4⋅x³ + 2⋅x²
    @assert mod(p1, p).lst == MutableLinkedList{Term128}(Term128(Int128(2),2), 
                                                            Term128(Int128(4),3), Term128(Int128(1),4), Term128(Int128(3),5), Term128(Int128(3),7))
    #mod(p1, 5) = 3⋅x⁷ + 2⋅x⁶ + x⁴ + 4⋅x³ + 2⋅x²
    @assert mod(p1, q).lst == MutableLinkedList{Term128}(Term128(Int128(2),2), 
                                                            Term128(Int128(4),3), Term128(Int128(1),4), Term128(Int128(2),6), Term128(Int128(3),7))
    #mod(p2, 7) = x⁷ + 2⋅x⁶ + 6⋅x⁵ + 6⋅x⁴ + 2⋅x³ + 6⋅x²
    @assert mod(p2, p).lst == MutableLinkedList{Term128}(Term128(Int128(6),2), Term128(Int128(2),3), 
                                                            Term128(Int128(6),4), Term128(Int128(6),5), Term128(Int128(2),6), Term128(Int128(1),7))
    #mod(p2, 5) = x⁷ + 2⋅x⁶ + x⁴ + 2⋅x³ + 3⋅x²
    @assert mod(p2, q).lst == MutableLinkedList{Term128}(Term128(Int128(3),2), 
                                                            Term128(Int128(2),3), Term128(Int128(1),4), Term128(Int128(2),6), Term128(Int128(1),7))
    #mod(p3, 7) = x⁷ + x⁶ + 2⋅x⁴
    @assert mod(p3, p).lst == MutableLinkedList{Term128}(Term128(Int128(2),4), Term128(Int128(1),6), Term128(Int128(1),7))
    #mod(p3, 5) = x⁷ + x⁶ + 2⋅x⁵
    @assert mod(p3, q).lst == MutableLinkedList{Term128}(Term128(Int128(2),5), Term128(Int128(1),6), Term128(Int128(1),7))

    println("mod_test_poly_sparse128 - PASSED")
end 

"""
Test capabilities of PolynomialSparse128 types at performing operations
on integers greater than 2⁶⁴ or less than -2⁶⁴. Note PolynomialSparse type cannot accept 
coefficients with values not in the range [-2⁶⁴, 2⁶⁴] and will "overflow".
"""
function bigint_test(u::Int128 = Int128(123456789012345678901234567890), l::Int128 = Int128(-123456789012345678901234567890))
    #Creat two polynomials with coefficients greater than 2⁶⁴ or less than -2⁶⁴
    p1 = PolynomialSparse128([Term128(u, 1), Term128(l, 5), Term128(l-100, 6)])
    p2 = PolynomialSparse128([Term128(l,1), Term128(u, 2), Term128(u, 3)])
    #add polynomials together (coefficients remain larger than what usual Sparse type can handle)
    p1 + p2
    #multiply polynomials togehter (some coefficients will remain larger or smaller than what usual Sparse type can handle)
    p1*p2
    #if no error is thrown, code has succesfully performed the operations where usual Sparse type would fail 
    @assert true 
    print("bigint_test - PASSED")
end 