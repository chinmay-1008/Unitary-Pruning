using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf
using LinearAlgebra


function run(; N = 6, a = 1, b = 1, w = 1, T = 10, dt = 1, thresh = 0)

    # ham = UnitaryPruning.central_spin(N, a, b, w)
    ket = KetBitString(N, 0) 

    k = floor(Int64(T/dt))

    generators, parameters = UnitaryPruning.central_spin_gen(N, a, b, w, k)

    # for ki in 1:k
    #     for (hi, coeff) in ham.ops
    #         push!(generators, Pauli(string(hi)))
    #         push!(parameters, coeff*dt)
    #     end
    # end
    # println("The Hamiltonian: ")
    # display(ham)
    # o = (1/sqrt(2))*(Pauli(N, X = [1]) + 1im*Pauli(N, Y = [1]))
    o_1 = PauliSum(N) + (1/sqrt(2))*(Pauli(N, X = [1]))
    o_2 = PauliSum(N) + (1im/sqrt(2))*(Pauli(N, Y = [1])) 
    o = o_1 + o_2
    println("------------------------------------------------------")

    nt = length(generators)
    println(nt, " ", k)
    ham_n = 3*(N-1)

    p = Int64(floor(k/4) + 1)*ham_n + 1

    println(ham_n, " ", p)
    insert!(generators, p, Pauli(N, X = [1]))
    insert!(parameters, p, π/2)
    println(generators)

    println("The observable: ")
    display(o)
    println("------------------------------------------------------")
 
    ei, nops, norm = UnitaryPruning.bfs_evolution_central(generators, parameters, o, ket, dt, T, thresh = thresh)
    @printf(" e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", real(ei), imag(ei), maximum(nops), norm, thresh)
 
    println("Expectation Value: ", ei)
    # println("Norm: ", norm)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println("Exact: ",(m[1]), " BFS: ", (ei), " Error: ",(m[1]) -  (ei) )
    return
end 

a = [0.176538, 0.174746, 0.169478, 0.161048, 0.149946]
b = [0.669628, 0.66283, 0.642846, 0.610871, 0.568759]
w = 2*π*3.98

run(N = 6, a = a, b = b, w = w, T = 2, dt = 0.5, thresh = -1)
# 