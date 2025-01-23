using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf


function run(; N=10, k=5,threshold = 1e-3, gamma = 0.5, lc = 1)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    # angles = [] 
    # e = [] 
    
    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = α,Jz = α, k=k)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_coeff_error(generators, parameters, PauliSum(o), ket, thresh = threshold, epsilon = epsilon)
    ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_evolution_diff(generators, parameters, PauliSum(o), ket, thresh = threshold, γ = gamma, lc = lc )

    # ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_angle_error(generators, parameters, PauliSum(o), ket, thresh = sigma, epsilon = epsilon)
 

    @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f γ: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, threshold, gamma)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei_n), " Error: ",real(m[1]) -  real(ei_n))

    return real(ei_n)
end


function run_temp()
    th = [-1, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3]
    for lc in 1:6
        for thre in th
            
            ep = 0.5e-3
            eps = [i*ep for i in 1:201]

            energy = []
            for ga in eps

                ei = run(N = 8, k = 10, threshold=thre, gamma = ga, lc = lc)
                push!(energy, ei)

            end

            # println(energy)
            # plot(eps, energy)
            # xlabel!("γ")
            # ylabel!("Expectation Value")
            # title!("(N=8; k=10; α=4π/32; thresh = 5e-3, l* = 3)")
            # savefig("test/ener_vs_gamma_N8_5e-3_3_new.pdf")
            writedlm("test/diss_data/exp_gam_N8_$thre_$lc.dat", energy)
        end
    end
    return
end

function run_samples()

    energies = []
    for time in 1:1000
        temp_energies = run_temp()

        println("Sample: ", time, ", completed")
        push!(energies, temp_energies)
        if time == 100
            writedlm("test/all_data/sam100_ev_N10_diss_1e-3_eps_1_10.dat", energies)
        end
    end

    writedlm("test/all_data/sam1000_ev_N10_diss_1e-3_eps_1_10.dat", energies)

    return
end


# @time run_samples()
run_temp()
# run(N = 10, k = 10, threshold = 1e-3, gamma = 1e-3, lc = 3)
# 