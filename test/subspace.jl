using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf


function run(; N=10, k=5,thresh = 1e-3)
   
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
    ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_subspace_expansion(generators, parameters, PauliSum(o), ket, thresh = thresh)
    # ei_n , nops_n, c_norm2_n = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh = thresh)



    @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, thresh)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei_n), " Error: ",real(m[1]) -  real(ei_n))

    return real(ei_n)
end

@time run(N = 3, k = 10, thresh = 1e-3)

