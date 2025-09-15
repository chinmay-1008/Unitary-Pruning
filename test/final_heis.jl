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
   
    ket = Ket(N, 0) 
    o = Pauli(N, Z=[1])

    # generators, parameters = UnitaryPruning.hzeisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx =1.0, Jy = 1.0,Jz = 1.0, k=k)
    # generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg_2D_open(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.fermi_hubbard_1D(o, t = 1, U = 6, k=k)
    generators, parameters, hammy = UnitaryPruning.fermi_hubbard_2D(o, t = 1, U = 2, k=k)

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution_weight(generators, parameters, PauliSum(o), ket, thresh=thresh, w_type = w_type, w = w)
    
    α = π/4

    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei), imag(ei), nops[end], c_norm2, thresh)

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


function run_ops(; N=10, k=5, thresh=1e-3, w_type = 0, w = 2)
   
    ket = Ket(N, 0) 
    o = Pauli(N, Z=[1])

    # generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg_2D_open(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.heisenberg(o, Jx =1.0, Jy = 1.0,Jz = 1.0, k=k)
    # generators, parameters = UnitaryPruning.heisenberg_2D(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    # generators, parameters = UnitaryPruning.fermi_hubbard_1D(o, t = 1, U = 6, k=k)
    generators, parameters = UnitaryPruning.fermi_hubbard_2D_new(o, t = 1, U = 2, k=k)

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution_weight(generators, parameters, PauliSum(o), ket, thresh=thresh, w_type = w_type, w = w)
    
    α = π/4

    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei), imag(ei), nops[end], c_norm2, thresh)

    return real(ei), nops[end], length(generators)
end

function run_weights_and_ops(run_weights_plot::Bool = true, run_ops_plot::Bool = true)
    L = 2
    N = 2 * L * L
    o = Pauli(N, Z=[1])
    new_set_k = [1, 2, 5, 10]
    # new_set_k = [1, 2, 3, 4, 5]

    # thresholds = [1e-4, 1e-3]
    thresholds = [-1, 1e-10, 1e-4, 1e-3]

    weights = [i for i in 1:2*N]

    for k in new_set_k
        println("k: ", k)
        # generators, parameters = UnitaryPruning.heisenberg(o, Jx=0.8, Jy=0.9, Jz=0.9, k=k)
        # generators, parameters = UnitaryPruning.heisenberg_2D_open(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
        generators, parameters = UnitaryPruning.fermi_hubbard_2D_new(o, t = 1, U = 2, k=k)

        # return
        # Placeholder for "exact" expectation value, set to 0
        # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
        # o_mat = Matrix(o)
        # m = diag(U'*o_mat*U)
        # display(m[1])
        m = [0.0]

        if run_weights_plot
            plt1 = plot(xlabel = "Weight Cutoff", ylabel = "Expectation Value",
                        title = "L = $L, N = $N, k=$k", grid = true, dpi = 300)#, legend = :bottomright)
        end

        if run_ops_plot
            plt2 = plot(xlabel = "Weight Cutoff", ylabel = "PP ops",
                        title = "L = $L, N = $N, k=$k", grid = true, dpi = 300)
        end

        for thresh in thresholds
            for w_type in (0, 1)  # 0: Pauli, 1: Majorana
                p_errors = Float64[]
                op_counts = Float64[]

                for w in weights
                    println("Weight Threshold: ", w)
                    # ev = run(N=N, k=k, thresh=thresh, w_type=w_type, w=w)

                    ev, nops, len_generators = run_ops(N=N, k=k, thresh=thresh, w_type=w_type, w=w)
                    push!(op_counts, nops)
                    push!(p_errors, abs((m[1]) - (ev)))

                end

                label = w_type == 0 ? "Pauli $thresh" : "Majorana $thresh"

                if run_weights_plot
                    plot!(plt1, weights, p_errors, label=label, lw=2, marker=:circle)
                end

                if run_ops_plot
                    plot!(plt2, weights, op_counts, label=label, lw=2, marker=:circle)
                end
            end
        end

        if run_weights_plot
            savefig(plt1, "test/hubbard_2d_weights_L$L-k$k.png")
        end
        if run_ops_plot
            savefig(plt2, "test/hubbard_2d_nops_L$L-k$k.png")
        end
    end
end


function eigenspectrum()
    N = 2
    N = 2*N*N
    k = 1
    o = Pauli(N, Z=[1])
    t = 1
    U = 12

    generators, parameters, hammy = UnitaryPruning.fermi_hubbard_2D(o, t = t, U = U, k = 1)
    e, v = eigen(Matrix(hammy))

    filename = "test/eigenspectrum_t$t-U$U.txt"

    open(filename, "w") do f
        for item in e
            println(f, item)
        end
    end

    return

end


function temp_run()
    N=5
    thresholds=[-1, 1e-6, 1e-5, 1e-4, 1e-3]
    ket = Ket(N, 0) 
    o = Pauli(N, Z=[1])

    k_list = [i for i in 1:20]
    plt = plot(xlabel = "k", ylabel = "# operators", dpi = 200,  title="Operator Growth vs k for N=$N")

    for thresh in thresholds
        nops_list = []
        println("Threshold: ", thresh)
        for k in k_list
            println("k: ", k)
            
            generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
            ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)       
            push!(nops_list, nops[end])
            # @printf(" e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", real(ei), imag(ei), nops[end], c_norm2, thresh)
        end
        plot!(plt, k_list, nops_list, label = "$thresh", lw=2, marker=:circle)
        display(nops_list)
    end
    savefig("test/nops_N$N.png")
    return 
end

function temp_ham()
    N = 3
    N = 2*N*N
    k = 8
    o = Pauli(N, Z=[1])
    t = 1
    U = 2
    ket = Ket(N, 0) 
    thresh = 1e-3
    w_type = 1
    w = 6

    generators, parameters = UnitaryPruning.fermi_hubbard_2D_new(o, t = t, U = U, k = k)
    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # display(m[1])

    # ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution_weight(generators, parameters, PauliSum(o), ket, thresh=thresh, w_type = w_type, w = w)
    
    if w_type == 1 

        println("Majorana Truncation")
    else
        println("Pauli Truncation")

    end
    @printf(" e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f weight: %3.1f\n", real(ei), imag(ei), nops[end], c_norm2, thresh, w)

    # for i in eachindex(generators)
    #     println(parameters[i])
    #     display(generators[i])
    # end
    # temp = KetSum(N)
end
# temp_run()

# temp_ham()
# run_weights()
# eigenspectrum()    
# run_num_ops()
run_weights_and_ops()