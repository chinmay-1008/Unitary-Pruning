using Distributed
@everywhere begin
    using UnitaryPruning
    using Plots
    using Statistics
    using Printf
    using Random
    using LinearAlgebra
    using SharedArrays
    using PauliOperators
end

function run(; N=10, k=5, thresh=1e-3, w = 2)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[5])

    # generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx =1.0, Jy = 1.0,Jz = 1.0, k=k)
    generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution_heisenberg(generators, parameters, PauliSum(o), ket, thresh=thresh, w = w)
    
    α = π/4

    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # error = real(m[1]) -  real(ei) 
    # println(real(m[1]), " ", real(ei), " Error: ", error)
    return real(ei)
end

function run_weights()
    N = 3
    k = 10
    N *= N
    o = Pauli(N, Z=[5])
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx =1.0, Jy = 1.0,Jz = 1.0, k=k)
    generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)


    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    o_mat = Matrix(o)
    m = diag(U'*o_mat*U)
    println("Exact Expectaition Value: ", m[1])

    errors = []
    weights = [i for i in 1:2*N]

    for i in weights
        println("Weight Threshold: ", i)
        e_v = run(N = N, k = k, thresh = -1, w = i)
        push!(errors, abs(real(m[1]) - real(e_v)))
    end
    display(errors)
    plot(weights, errors;
        label = "Error",
        xlabel = "Majorana Weight Cutoff",
        ylabel = "Abs Error",
        title = "Erros vs Majorana Weight for N = $N, k=$k",
        lw = 2,
        marker = :circle,
        legend = :topright,
        grid = true)
    savefig("test/fin_heis2d_weights_N$N-k$k-majorana_2.pdf")

end

run_weights()
# N = 3
# N = N^2
# k = 1
# o = Pauli(N, Z=[1])

# generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
# # display(generators)

# L = Int(sqrt(N))
# index(i, j) = (mod1(i, L) - 1) * L + mod1(j, L)


