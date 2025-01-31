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

function msd(N, cj)
    total_time = length(cj)
    msd_val = []
    for time in 1:total_time
        temp_c = Vector{Float64}(cj[time])
        msd_dt = 0
        msd_1 = 0
        msd_2 = 0
        for j in 1:N
            msd_1 += temp_c[j]*(j^2)
            msd_2 += temp_c[j]*j
        end
        msd_dt = msd_1 - (msd_2^2)
        push!(msd_val, msd_dt)
        # println(temp_c)
    end
    
    return msd_val
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
    c_j = UnitaryPruning.bfs_evolution_new_diff(generators, parameters, o, hj, k, thresh =0, γ = γ, lc = lc)

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
    val = msd(N, c_j)
    time_steps = [dt*j for j in 1:k]
    plot(time_steps, val)
    xlabel!("Time")
    ylabel!("MSD")
    title!("(N=9; dt=0.25; thresh = 0; γ=0.03, l* = 2)")
    savefig("test/new_diff_N9_paper_0.pdf")
end    


run(N = 9, Jx = 1.4, Jz = 0.9045, dt = 0.25, T = 10, γ = 0.03, lc = 2)
