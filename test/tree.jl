using UnitaryPruning
using PauliOperators
using LinearAlgebra
using DelimitedFiles

function run(; N=10, k=5,threshold = 1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])
    # o = Pauli(N, X=[2], Y=[1], Z = [3])

    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8/2, Jy = 0.9/2,Jz = 0.9/2, k=k)
    # ei_n , nops_n, c_norm2_n, dict = UnitaryPruning.tree_evolution(generators, parameters, PauliSum(o), ket, thresh = threshold)
    otoc, dict = UnitaryPruning.tree_evolution_otoc(generators, parameters, PauliSum(o), ket, thresh = threshold)
 

    # @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, threshold)

    println(otoc)
    evol = []
    for data in dict
        temp_pauli = []
        temp_coeff = []
        for (oi, coeff) in data.ops
            push!(temp_pauli, string(oi))
            push!(temp_coeff, coeff)
        end
        push!(evol, temp_pauli)
        push!(evol, temp_coeff)
    end
    time_steps = [ i for i in 1:k]
    plot(time_steps, otoc)
    xlabel!("Time")
    ylabel!("OTOC")
    title!("N=$N; thresh = $threshold")
    savefig("test/tree_otoc$N-$threshold.pdf")
 
    # writedlm("test/tree_evolN$N-$threshold.dat", evol)
    return
end

@time run(N = 10, k = 20, threshold = 1e-3)
# a = Pauli(6)
# display(a)
# println(string(a))