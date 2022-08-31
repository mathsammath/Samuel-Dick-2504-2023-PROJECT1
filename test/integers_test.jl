#############################################################################
#############################################################################
#
# This file contains defines unit tests for integer operations
#                                                                               
#############################################################################
#############################################################################


"""
Tests Euclidean GCD algorithm for integers.
"""
function test_euclid_ints(;N::Int = 10^4)
    Random.seed!(0)
    for _ in 1:N
        n1 = rand(1:10^6)
        n2 = rand(1:10^6)
        g = euclid_alg(n1,n2)
        @assert mod(n1,g) == 0 &&  mod(n2,g) == 0
    end
    println("test_euclid_ints - PASSED")
end

"""
Tests the extended Euclidean GCD algorithm for integers.
"""
function test_ext_euclid_ints(;N::Int = 10^4)
    Random.seed!(0)
    for _ in 1:N
        n1 = rand(1:10^6)
        n2 = rand(1:10^6)
        g, s, t = ext_euclid_alg(n1,n2)
        @assert g == s*n1 + t*n2
        @assert mod(n1,g) == 0 &&  mod(n2,g) == 0
    end
    println("test_ext_euclid_ints - PASSED")
end

"""
Tests the computation of inverse mod for integers.
"""
function test_inverse_mod_ints(;prime::Int=101,N::Int=10^4)
    Random.seed!(0)
    for _ in 1:N
        n = rand(1:10^6)
        if mod(n, prime) == 0
            continue
        end
        im = int_inverse_mod(n,prime) 
        @assert mod((n*im), prime) == 1
    end
    println("test_inverse_mod_ints - PASSED")
end