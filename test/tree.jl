using UnitaryPruning
using PauliOperators
using LinearAlgebra
using DelimitedFiles

function run(; N=10, k=5,threshold = 1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])
    # o = Pauli(N, X=[2], Y=[1], Z = [3])

    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    ei_n , nops_n, c_norm2_n, dict = UnitaryPruning.tree_evolution(generators, parameters, PauliSum(o), ket, thresh = threshold)
 

    @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, threshold)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei_n), " Error: ",real(m[1]) -  real(ei_n))
    # display(dict)
    # pauli_evol = []
    evol = []
    for data in dict
        temp_pauli = []
        temp_coeff = []
        for (oi, coeff) in data.ops
            push!(temp_pauli, string(oi))
            push!(temp_coeff, coeff)
        end
        push!(evol, temp_pauli)
        push!(evol, temp_coeff)
    end
    writedlm("test/tree_evolN$N-$threshold.dat", evol)
    return
    # println(dict[1])
end

run(N = 5, k = 5, threshold = 0)
# a = Pauli(6)
# display(a)
# println(string(a))