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

function run(; N=10, k=5, thresh=1e-3, w_type = 0, w = 2)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    # generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx =1.0, Jy = 1.0,Jz = 1.0, k=k)
    # generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    generators, parameters = UnitaryPruning.fermi_hubbard_1D(o, t = 1, U = 8, k=k)

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution_weight(generators, parameters, PauliSum(o), ket, thresh=thresh, w_type = w_type, w = w)
    
    α = π/4

    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)

    # U = UnitaryPruning.build_time_evolution_matrix_fast(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # error = real(m[1]) -  real(ei) 
    # println(real(m[1]), " ", real(ei), " Error: ", error)


    # ψ = zeros(ComplexF64, 2^N); ψ[1] = 1.0
    # ψp = U * ψ
    # e = dot(ψp, Matrix(o) * ψp)
    # error = real(e) -  real(ei) 
    # println(real(e), " ", real(ei), " Error: ", error)

    return real(ei)
end

function run_weights()
    L = 5

    k = 10
    N = 2 * L 
    o = Pauli(N, Z=[1])
    thresh = 1e-4
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx =1.0, Jy = 1.0,Jz = 1.0, k=k)
    # generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    generators, parameters = UnitaryPruning.fermi_hubbard_1D(o, t = 1 , U = 8, k=k)


    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    m = [0]
    println("Exact Expectaition Value: ", m[1])

    weights = [i for i in 1:2*N]

    w_type = 0
    p_errors = []

    for i in weights
        println("Weight Threshold: ", i)
        e_v = run(N = N, k = k, thresh = thresh, w_type = w_type, w = i)
        push!(p_errors, abs(real(m[1]) - real(e_v)))
        # push!(p_errors, real(e_v))

    end
    display(p_errors)
    # type = ""
    # if w_type == 0
    #     type = "Pauli"
    # elseif w_type == 1
    #     type = "Majorana"
    # end

    # plot(weights, errors;
    #     label = "$type Error",
    #     xlabel = "$type Weight Cutoff",
    #     ylabel = "Abs Error",
    #     title = "Erros vs $type Weight for L = $L, N = $N, k=$k",
    #     lw = 2,
    #     marker = :circle,
    #     legend = :topright,
    #     grid = true)

    
    w_type = 1
    m_errors = []

    for i in weights
        println("Weight Threshold: ", i)
        e_v = run(N = N, k = k, thresh = thresh, w_type = w_type, w = i)
        push!(m_errors, abs(real(m[1]) - real(e_v)))
        # push!(m_errors, real(e_v))

    end
    # display(errors)
    # type = ""
    # if w_type == 0
    #     type = "Pauli"
    # elseif w_type == 1
    #     type = "Majorana"
    # end

    plot(weights, [p_errors, m_errors],
        label = ["Pauli_Exp" "Majorana_Exp"],
        xlabel = "Weight Cutoff",
        ylabel = "Expectation Value",
        title = "L = $L, N = $N, k=$k, thresh = $thresh",
        lw = 2,
        marker = :circle,
        legend = :topright,
        grid = true)

    savefig("test/hubbard_1d_weights_L$L-k$k-th$thresh-2.pdf")

end

run_weights()
# @time run(N = 10, k = 10, thresh = -1, w_type = 0, w=8)
# N = 2
# N = 2*N
# k = 1
# o = Pauli(N, Z=[1])

# generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
# # display(generators)
# generators, parameters = UnitaryPruning.fermi_hubbard_1D(o, t = 0.1, U = 1, k = 1)

# display(hammy)
# p = UnitaryPruning.jw_transform(o, 1)
# display(p' * p)
# L = Int(sqrt(N))
# index(i, j) = (mod1(i, L) - 1) * L + mod1(j, L)


