using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf


function run(; N=10, k=5,sigma = 1e-3, layer = false)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    angles = [] 
    e = []
    e_mat = [] 
    n = 30
    for i in 1:2*n
        # i = 4
        α = i / n 
        generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = 1.0,Jz = 1.0, k=k)
        # ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
        ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh = sigma)
        

        @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f sigma: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, sigma)

        # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
        # o_mat = Matrix(o)
        # m = diag(U'*o_mat*U)
        # println(real(m[1]), " ", real(ei_n), " Error: ",real(m[1]) -  real(ei_n))
        push!(angles, α)
        push!(e, real(ei_n))
        # push!(e_mat, real(m[1]))
    end
 
    # println(angles)
    # println()
    # println(e)
    # println()
    # println(e_mat)
    scatter(angles, e)
    # scatter(angles, [e_mat, e], label = ["Exact" "BFS"])
    # p2 = scatter(angles, e, label = "BFS")

    # # plot(p1, p2, layout=(1, 2))
    xlabel!("J_x")
    ylabel!("Expectation Value")
    title!("E{Z_1} vs J_x (N=8; k=10; J_y = 1.0; J_z = 1.0)")
    savefig("test/ener_vs_Jx_8_f_det.pdf")
    println(e)
    return e
end


function run_exact(; N=10, k=5,sigma = 1e-3, layer = false)
   
    o = Pauli(N, Z=[1])
    angles = [] 
    e_mat = [] 
    n = 30
    for i in 1:2*n
        α = i / n 
        generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = 1.0,Jz = 1.0, k=k)
        
        U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
        o_mat = Matrix(o)
        m = diag(U'*o_mat*U)

        push!(angles, α)
        push!(e_mat, real(m[1]))
    end
    return [e_mat, angles]
end

function run_samples()
    println("Exact Calculating")
    exact = run_exact(N = 8, k = 10, sigma = 1e-3, layer = false)
    writedlm("test/sam1000_ev_N8_Jx_exact.dat", exact)
    println("Exact data calculated")

    energies = []
    for time in 1:1000
        temp_energies = run(N = 8, k = 10, sigma = 1e-3, layer = false)

        println("Sample: ", time, ", completed")
        push!(energies, temp_energies)
    end

    writedlm("test/sam1000_ev_N8_Jx.dat", energies)

    return
end

@time run(N = 8, k = 10, sigma = 1e-2, layer = false)
# @time run_samples()