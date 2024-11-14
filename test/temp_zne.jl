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
    
    # for i in [(i-1)*2 for i in 1:9]
    # for i in 9:9
    i = 9
    α = i * π / 32 
    generators, parameters = UnitaryPruning.get_unitary_sequence_1D(o, α=α, k=k)
    
    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
    
    push!(e, ei)
    push!(angles, α)
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f\n", α, real(ei), imag(ei), maximum(nops), c_norm2)
    # end
    
    # plot(angles, real(e), marker = :circle)
    # xlabel!("Angle")
#    ylabel!("expectation value")
#    title!{X_{13,29,31}, Y_{9,30}, Z_{8,12,17,28,32}}
    # savefig("plot_1d_n6_bfs.pdf")
    return real(ei)
end

thresholds = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4]
# thresholds = [0]
ener = []
for th in thresholds
    # println(th)
    en = run(k=10, N=12, thresh=th)
    push!(ener, en)
end
println(ener)
# plot(thresholds, ener, marker = :circle)
# savefig("test/ener_vs_thre.png")