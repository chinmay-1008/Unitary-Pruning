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
        generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = 0.8,Jz = 0.6, k=k)
        ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
        

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
    title!("E{Z_1} vs J_x (N=127; k=10; J_y = 0.8; J_z = 0.6)")
    savefig("test/ener_vs_Jx_127_f.pdf")
    return
end


@time run(N = 127, k = 10, sigma = 1e-3, layer = false)