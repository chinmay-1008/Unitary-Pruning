using UnitaryPruning
using PauliOperators
using Plots
using DelimitedFiles
using Printf


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

function generator_ising(o::Pauli{N}; Jx, Jz, k) where N 
    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])
    # Loop over sites
    for ki in 1:k 
        for i in 1:N
            push!(generators, Pauli(N, X=[i]))
            push!(parameters, Jx)
            push!(generators, Pauli(N, Z=[i]))
            push!(parameters, Jz)
            push!(generators, (Pauli(N, Z=[(i-1+N-1)%N + 1, i])))
            push!(0.5)
            push!(generators, Pauli(N, Z=[i, (i%N) + 1]))
            push!(0.5)
        end
    end

    return generators, parameters
end

function run()

end    

