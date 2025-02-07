using UnitaryPruning
using PauliOperators
using LinearAlgebra
using DelimitedFiles
using Plots

function run(; N=10,threshold = 1e-3, dt=0.1, T=10)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])
    # o = Pauli(N, X=[2], Y=[1], Z = [3])
    k = Int64(floor(T/dt))
    x = 0.8
    y = 0.9
    z = 0.9
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = x*dt, Jy = y*dt,Jz = z*dt, k=k)
    # ei_n , nops_n, c_norm2_n, dict = UnitaryPruning.tree_evolution(generators, parameters, PauliSum(o), ket, thresh = threshold)
    otoc, dict = UnitaryPruning.tree_evolution_otoc(generators, parameters, PauliSum(o), ket, thresh = threshold)
 

    # @printf("α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, threshold)

    # println(otoc)
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
    time_steps = [ dt*i for i in 1:k+1]

    plot(time_steps, otoc, titlefontsize = 5)
    xlabel!("Time")
    ylabel!("OTOC")
    title!("N=$N; Jx = $x, Jy = $y, Jz = $z, thresh = $threshold; dt = $dt")
    savefig("test/tree_otoc$N-$threshold-$dt-$T-new-abc-2.pdf")
 
    # writedlm("test/tree_evolN$N-$threshold.dat", evol)
    return
end

@time run(N = 5, threshold = 1e-3, dt = 0.1, T = 100)
# a = Pauli(6)
# display(a)
# println(string(a))