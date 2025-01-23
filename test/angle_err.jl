using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf


function run(; N=10, k=5,sigma = 1e-3,epsilon = 1e-3, layer = false)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    # angles = [] 
    # e = [] 
    
    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = α,Jz = α, k=k)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.flat_prob_kill_bfs(generators, parameters, PauliSum(o), ket, threshold = sigma)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs_flat(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs_renormalized(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs_coeff_error(generators, parameters, PauliSum(o), ket, sigma = sigma, epsilon = epsilon,layer = layer)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_coeff_error(generators, parameters, PauliSum(o), ket, thresh = sigma, epsilon = epsilon)
    ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_angle_error(generators, parameters, PauliSum(o), ket, thresh = sigma, epsilon = epsilon)
 

    # @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f sigma: %3.10f epsilon: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, sigma, epsilon)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei_n), " Error: ",real(m[1]) -  real(ei_n))

    return real(ei_n)
end


function run_temp()
    ep = 1e-3
    eps = [i*ep for i in 1:1000]
    # temp_std = [i*st for i in 30:10:100]
    # for j in temp_std
    #     push!(stds, j)
    # end
    # println(stds)

    energy = []
    for epsilon in eps

        ei = run(N = 5, k = 10, sigma=0, epsilon=epsilon, layer = false)
        push!(energy, ei)

    end

    # println(energy)
    # scatter(stds, energy, smooth=true)
    # xlabel!("σ")
    # ylabel!("Expectation Value")
    # title!("E{Z_1} vs σ (N=10; k=10; α=4π/32)")
    # savefig("test/ener_vs_sigma_true.pdf")
    return energy
end

function run_samples()

    energies = []
    for time in 1:1000
        temp_energies = run_temp()

        println("Sample: ", time, ", completed")
        push!(energies, temp_energies)
    end

    writedlm("test/all_data/sam100_ev_f_N5_coeff_an_large.dat", energies)

    return
end


@time run_samples()
# run(N = 5, k = 10, sigma = 0, epsilon = 1e-1, layer = false)