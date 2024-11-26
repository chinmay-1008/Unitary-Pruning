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

function run(; N=10, k=5, thresh=1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    angles = [] 
    e = [] 
    
    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = α,Jz = α, k=k)
    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
    
    push!(e, ei)
    push!(angles, α)
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei), " Error: ",real(m[1]) -  real(ei) )
    return real(ei)
end


function run_local_folding()

    N = 10
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])
    k = 10
    i = 4
    α = i * π / 32 
    thresh=1e-3
    scaling_factor = [j for j in 0:20]
    energy = []

    for additional in scaling_factor

        # generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
        generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = α,Jz = α, k=k)

        generators_lf, parameters_lf = UnitaryPruning.local_folding(generators, parameters, additional)
        # println(generators_lf)
        println(length(generators)," " ,length(generators_lf))

        ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
        @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.5f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)
        # println(nops)
        ei_lf , nops_lf, c_norm2_lf = UnitaryPruning.bfs_evolution(generators_lf, parameters_lf, PauliSum(o), ket, thresh=thresh)

        println("Additional Pauli Terms: ", 2 * additional)
        @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei_lf), imag(ei_lf), maximum(nops_lf), c_norm2_lf, thresh)
        # println(nops_lf)
        push!(energy, real(ei_lf))

    end
    println(energy)
    scatter(scaling_factor, energy, smooth=true)
    xlabel!("λ")
    ylabel!("Expectation Value")
    title!("E{Z_1} vs # UU'(N=10; k=10; α = 4π/32)")
    savefig("test/ener_vs_thre_hei_lf.pdf")

    return 
end

function run_threshold()
    # Getting expectation values at different values of k 
  
    # thresholds = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4]
    # thresholds = [9e-5, 8e-5, 6e-5, 5e-5, 4e-5, 3e-5, 2e-5, ]
    # thresholds = [5.0e-4, 5.1e-4, 5.2e-4, 5.3e-4, 5.5e-4, 5.6e-4, 5.7e-4, 5.8e-4, 5.9e-4, 6.0e-4, 6.1e-4, 6.2e-4, 6.3e-4, 6.4e-4, 6.5e-4, 6.6e-4, 6.7e-4, 6.8e-4, 6.9e-4, 7.0e-4]
    thresholds = [i*(1e-4) for i in 1:20]
    ener = []
    for th in thresholds
        # println(th)
        en = run(N=10, k=10, thresh=th)
        push!(ener, en)
        # println(th, " ", en)
    end

    println(ener)
    scatter(thresholds, ener, smooth=true)
    xlabel!("Threshold")
    ylabel!("Expectation Value")
    title!("E{Z_1} vs threshold (N=10; k=10; α=4π/32)")
    savefig("test/ener_vs_thre_hei.pdf")
    return
end


# run(N = 10, k = 1, thresh = 1e-3)

run_local_folding()

# run_threshold()
