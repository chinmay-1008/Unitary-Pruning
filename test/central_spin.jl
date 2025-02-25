using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf
using LinearAlgebra


function run(; N = 6, a = 1, b = 1, w = 1, k = 10, T = 1, thresh = 0)

    # ham = UnitaryPruning.central_spin(N, a, b, w)
    ket_0 = KetBitString(N, 0)
    ket_1 = KetBitString(N, 1)

    new_state = SparseKetBasis(N)

    new_state[ket_0] = 1/sqrt(2)
    new_state[ket_1] = 1/sqrt(2)
    
    # display(new_state)

    # o_1 = (1/sqrt(2))*(Pauli(N, X = [1]))
    # o_2 = (1/sqrt(2))*(Pauli(N, Y = [1])) 
    # o = o_1 + o_2
    
    # ev_state = o * new_state
    # display(ev_state)

    # function dot_pdt(bra_state, ket_state)
    #     val = 0
    #     for (ket, coeff_ket) in ket_state
    #         for (bra, coeff_bra) in bra_state
    #             if ket.v == bra.v
    #                 val += coeff_bra' * coeff_ket
    #             end
    #         end
    #     end
    #     return val
    # end 
    # a = dot_pdt(new_state, ev_state)
    # println(a)
    # println(expectation_value(o, new_state))

    # return
    # dt = 0.01

    # k = Int64(T/dt)
    # T = floor(Int64(k*dt))
    dt = T/k
    println("dt: ", dt)
    generators, parameters = UnitaryPruning.central_spin_gen(N, a, b, w, k)

    o_1 = (1/sqrt(2))*(Pauli(N, X = [1]))
    o_2 = (1/sqrt(2))*1im*(Pauli(N, Y = [1])) 
    o = o_1 + o_2

    # o = PauliSum(N) + Pauli(N, Z = [1])
    # println("------------------------------------------------------")

    nt = length(generators)
    println("Num of generators: ", nt, " Trotter Steps: ", k)
    ham_n = 3*(N-1)

    # p_1 = Int64(floor(k/4))*ham_n + 1

    # println(ham_n, " ", p_1)
    # insert!(generators, p_1, Pauli(N, X = [1]))
    # insert!(parameters, p_1, π/2)

    # p_2 = (Int64(floor((3*k)/4)))*ham_n + 2

    # println(ham_n, " ", p_2)
    # insert!(generators, p_2, Pauli(N, X = [1]))
    # insert!(parameters, p_2, π/2)

    # pi_pulse = [p_1, p_2]
    pi_pulse = []

    # println(generators)
    println("The observable: ")
    display(o)
    println("------------------------------------------------------")
    println("Beginning the evolution: ")
    ei, nops, norm, ei_time, dicts = UnitaryPruning.bfs_evolution_central(generators, parameters, o, new_state, dt, pi_pulse, thresh = thresh)
    @printf(" e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", real(ei), imag(ei), maximum(nops), norm, thresh)
 
    println("Expectation Value: ", ei)

    time_step = [i*dt for i in 1:k]
    plot(time_step, ei_time)
    xlabel!("Time")
    ylabel!("Abs(Exp Val)")


    if length(pi_pulse) == 0
        title!("N=$N, k=$k, dt=$dt; no π-pulse")
        savefig("test/temp-central_$N-2pi-$k-$T-skb-nopi.pdf")

    else 
        title!("N=$N, k=$k, dt=$dt; π-pulse at t = T/4, 3T/4")
        t1 = T/4
        t2 = 3T/4
        ylims = [minimum(ei_time), maximum(ei_time)]  # Adjust y-limits for visibility
    
        # Overlay the vertical dashed red lines
        plot!([t1, t1], ylims, color=:red, linestyle=:dash, linewidth=1, label="T/4")
        plot!([t2, t2], ylims, color=:red, linestyle=:dash, linewidth=1, label="3T/4")
    
        savefig("test/temp-central_$N-2pi-$k-$T-skb-pi.pdf")
    end



    # evol = []
    # for data in dicts
    #     temp_pauli = []
    #     temp_coeff = []
    #     for (oi, coeff) in data.ops
    #         push!(temp_pauli, string(oi))
    #         push!(temp_coeff, coeff)
    #     end
    #     push!(evol, temp_pauli)
    #     push!(evol, temp_coeff)
    # end
    # writedlm("test/central_evolN$N-$k-$T-skb-nopi.dat", evol)
    # println(dicts)
    # savefig("test/central_$N-2pi-$k-$T-nuc6.pdf")
    # display(ei_time)
    # c_t = ei_time.*ei_time
    # one_tangling = []
    # for i in c_t
    #     plus = 1 + i
    #     minus = 1 - i

    #     temp = -(plus/2)*log2(plus/2) - (minus/2)log2(minus/2)
    #     push!(one_tangling, temp)
    # end

    # plot(time_step, one_tangling)
    # xlabel!("Time")
    # ylabel!("One-tangling power")
    # title!("N=$N, k=$k, dt=$dt; π-pulse at t=T/4, 3T/4")

    # println("Norm: ", norm)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println("Exact: ",(m[1]), " BFS: ", (ei), " Error: ",(m[1]) -  (ei) )
    return
end 

a = [0.176538, 0.174746, 0.169478, 0.161048, 0.149946]
b = [0.669628, 0.66283, 0.642846, 0.610871, 0.568759]
w = 2*π*3.98

# run(N = 6, a = a, b = b, w = w, k = 500, T = 10, thresh = -1)

run(N = 2, a = [1], b = [1], w = 1, k = 5, T = 1, thresh = -1)
