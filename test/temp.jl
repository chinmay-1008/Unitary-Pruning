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
sin = PauliSum(5)
l = ReentrantLock()

@Threads.threads for (i, j) in collect(temp_pauli.ops)
    cos = PauliSum(5) + h[1]

    for m in 1:4
        tan = PauliSum(5)
        for (a, b) in cos.ops
            sum!(tan, a*b*0.5)
        end

        sum!(cos, tan)

    end
    lock(l) do 
        sum!(sin, cos)
    end

end

display(sin)


