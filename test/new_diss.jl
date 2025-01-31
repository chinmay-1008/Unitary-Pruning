using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf
using LinearAlgebra


function tilted_ising(N, Jx, Jz)
    H = PauliSum(N)
    hj = []
    for i in 1:N
        temp = PauliSum(N)

        temp += Jx * Pauli(N, X=[i])
        temp += Jz * Pauli(N, Z=[i])
        temp += 0.5 * (Pauli(N, Z=[(i-1+N-1)%N + 1, i]) + Pauli(N, Z=[i, (i%N) + 1]))
        H += temp
        push!(hj, temp)
    end 
    return H, hj
end 

function generator_ising(N; Jx, Jz, dt, k)
    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])
    # count = 0
    # Loop over sites
    for ki in 1:k 
        for i in 1:N
            push!(generators, Pauli(N, X=[i]))
            push!(parameters, Jx*dt)
            push!(generators, Pauli(N, Z=[i]))
            push!(parameters, Jz*dt)
            push!(generators, (Pauli(N, Z=[(i-1+N-1)%N + 1, i])))
            push!(parameters, 0.5*dt)
            push!(generators, Pauli(N, Z=[i, (i%N) + 1]))
            push!(parameters, 0.5*dt)
        end
    end

    return generators, parameters
end



function run(;N = 10, Jx = 1.0, Jz = 1.0, dt = 0.1, T = 10, γ = 0, lc = 0)
    H_i, hj = tilted_ising(N, Jx, Jz)
    o = hj[Int64((N+1)/2)] 

    println("Hamiltonian: ")
    display(H_i)
    println(" ")
    println("State: ")
    display(o)

    k = T/dt
    generators, parameters = generator_ising(N, Jx = Jx, Jz = Jz, dt = dt, k = 1)
    c_j = UnitaryPruning.bfs_evolution_new_diff(generators, parameters, o, hj, k, thresh = 0, γ = γ, lc = lc)

    # o_transformed = PauliSum(N)

    # for (oi, coeff) in o.ops
    #     temp_oi = PauliSum(N)
    #     println(" ")
    #     println("Initial Operator: ")
    #     display(oi)
    #     oi += temp_oi
    #     temp_o, n_ops = UnitaryPruning.bfs_evolution_new_diff(generators, parameters, oi, hj, coeff, k, thresh = 0, γ = γ, lc = lc)
    #     println("Number of operators after evolution: ", maximum(n_ops))
    #     o_transformed += temp_o
    # end 

    println(" ")
    println(c_j)
    # println("Evolved State: ")
    # display(o_transformed)
    # println("Max operators: ", max(n_ops))
    # c = adjoint(hj[1]) * o_transformed
    # println(length(o_transformed), " ", length(c)," ", length(hj[1]))
    # for (oi, coeff) in o
    # final_o, final_nops = UnitaryPruning.bfs_evolution_new_diff(generators, parameters,)
    # o_diss = dissipation(o_transformed, γ = 0.01, lc = 2)
    # display(o_diss)
end    


run(N = 5, Jx = 1.4, Jz = 0.9045, dt = 0.25, T = 10, γ = 0.025, lc = 2)
