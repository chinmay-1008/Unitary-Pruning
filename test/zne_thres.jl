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

function run(; N=6, k=10, thresh=1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    angles = [] 
    e = [] 
    
    i = 8
    α = i * π / 32 
    generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
    
    push!(e, ei)
    push!(angles, α)
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.5f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)

    return real(ei)
end

function run_zne()
    # Getting expectation values at different values of k 

    # thresholds = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4]
    thresholds = [9e-5, 8e-5, 6e-5, 5e-5, 4e-5, 3e-5, 2e-5, ]
    # thresholds = [5.0e-4, 5.1e-4, 5.2e-4, 5.3e-4, 5.5e-4, 5.6e-4, 5.7e-4, 5.8e-4, 5.9e-4, 6.0e-4, 6.1e-4, 6.2e-4, 6.3e-4, 6.4e-4, 6.5e-4, 6.6e-4, 6.7e-4, 6.8e-4, 6.9e-4, 7.0e-4]
    ener = []
    for th in thresholds
        # println(th)
        en = run(k=15, N=127, thresh=th)
        push!(ener, en)
        println(th, " ", en)
    end
    println(ener)
    # plot(thresholds, ener, marker = :circle)
    # xlabel!("Threshold")
    # ylabel!("Expectation Value")
    # title!("Expectation Value trend with threshold N = 127, k = 15")
    # savefig("test/ener_vs_thre_init_127.png")
    return
end

function local_folding(generators::Vector{Pauli{N}}, angles, scaling_factor) where {N}
    # Adding U'U at rondom places in the generators

    places = rand(1:(length(generators)-1), scaling_factor)
    places = sort(places, rev = true)

    temp_gens = deepcopy(generators)
    temp_angles = deepcopy(angles)

    for j in places
        add_pauli = temp_gens[j]
        add_angle = temp_angles[j]
        # temp_gens = vcat(temp_gens[1:j-1], add_pauli, add_pauli, temp_gens[j:end])
        # temp_angles = vcat(temp_angles[1:j-1], add_angle, -add_angle, temp_angles[j:end])

        insert!(temp_gens, j+1, -add_pauli)
        insert!(temp_gens, j+2, add_pauli)

        insert!(temp_angles, j+1, add_angle)
        insert!(temp_angles, j+2, add_angle)
      

    end
    return temp_gens, temp_angles
end

function run_lf()

    N = 15
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])
    k = 10
    i = 8
    α = i * π / 32 
    thresh=1e-3
    additional = 100
    generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
    generators_lf, parameters_lf = local_folding(generators, parameters, additional)

    println(length(generators)," " ,length(generators_lf)," ", length(parameters)," ", length(parameters_lf))

    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.5f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)
    # println(nops)
    ei_lf , nops_lf, c_norm2_lf = UnitaryPruning.bfs_evolution(generators_lf, parameters_lf, PauliSum(o), ket, thresh=thresh)

    println("Additional Pauli Terms: ", 2 * additional)
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.5f\n", α, real(ei_lf), imag(ei_lf), maximum(nops_lf), c_norm2_lf, thresh)
    # println(nops_lf)

    return
end

run_lf()
# run_zne()