using Random
using Distributions


function tilted_ising(N, Jx, Jz)
    H = PauliSum(N)
    hj = []
    for i in 1:N
        temp = PauliSum(N)

        temp += Jx * Pauli(N, X=[i])
        temp += Jz * Pauli(N, Z=[i])
        temp += 0.5 * (Pauli(N, Z=[(i-1+N-1)%N + 1, i]) + Pauli(N, Z=[i, (i%N) + 1]))
        H += temp
        push!(hj, temp)
    end 
    return H, hj
end 

temp_pauli, h = tilted_ising(5, 1, 1)

temp_pauli_1, h_1 = tilted_ising(5, 1.32, 2.65)

# display(temp_pauli)
sin_b = PauliSum(5)
l = ReentrantLock()
cos_b = PauliSum(5)
@Threads.threads for (i, j) in collect(temp_pauli.ops)
    cos_b = PauliSum(5) + h[1]

    for m in 1:4
        tan_b = PauliSum(5)
        for (a, b) in cos_b.ops
            sum!(tan_b, a*b*0.5)
        end

        sum!(cos_b, tan_b)

    end
    lock(l) do 
        sum!(sin_b, cos_b)
    end

end

display(sin_b)


