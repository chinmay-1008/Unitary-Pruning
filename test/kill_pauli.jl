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
    
    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = α,Jz = α, k=k)
    ei_n , nops_n, c_norm2_n = UnitaryPruning.stochastic_bfs(generators, parameters, PauliSum(o), ket, sigma = sigma, layer = layer)
    

    # @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f sigma: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, sigma)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei_n), " Error: ",real(m[1]) -  real(ei_n))

    return real(ei_n)
end


function run_temp()
    stds = [i*1e-3 for i in 1:20]
    energy = []
    for sigma in stds

        ei = run(N = 50, k = 10, sigma=sigma, layer = false)
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
    for time in 1:100
        temp_energies = run_temp()

        println("Sample: ", time, ", completed")
        push!(energies, temp_energies)
    end

    writedlm("test/sam_ev_layer_f_N50_.dat", energies)

    return
end


run_samples()
# run(N = 5, k = 1, sigma = 1e-3, layer = true)