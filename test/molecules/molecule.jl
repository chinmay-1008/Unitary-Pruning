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
    using JSON
end


function molecule()

    # For H2 molecue
    # pauli_ops = ["IIII", "IIIZ", "IIZI", "IZII", "ZIII", "IIZZ", "IZIZ", "ZIIZ",
    #             "YYYY", "XXYY", "YYXX", "XXXX", "IZZI", "ZIZI", "ZZII"]

    # coeffs=[-0.81217061,  0.17141283, -0.22343154,  0.17141283,
    #         -0.22343154,  0.12062523,  0.16868898,  0.16592785,
    #         0.04530262,  0.04530262,  0.04530262,  0.04530262,
    #         0.16592785,  0.17441288,  0.12062523]

    data = JSON.parsefile("test/molecules/h2_qubit_ham.json")
    pauli_ops = data["paulis"]
    coeffs = ComplexF64[]
    for (re, im) in zip(data["coeffs_real"], data["coeffs_imag"])
        push!(coeffs, complex(re, im))
    end

    N = length(pauli_ops[1])
    ket = KetBitString(N, 10) 
    display(ket)
    H = PauliSum(N)
    for i in eachindex(pauli_ops)
        sum!(H, Pauli(pauli_ops[i])*coeffs[i])
    end
    display(H)
    return
end

function ansatz()
    data = JSON.parsefile("test/molecules/parametric_generators.json")
    generators = []
    println("Ansatz")
    for g in data
        terms = PauliSum(length(g["terms"][1]["pauli"]))  # assumes all same length
        param = g["param"]
        for term in g["terms"]
            coeff = complex(term["real"], term["imag"])
            sum!(terms, Pauli(term["pauli"]) * coeff)
        end
        push!(generators, terms)
    end
    println(generators)
    return   # List of (param_name, PauliSum)
end

function run(; N=10, k=5, thresh=1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    angles = [] 
    e = [] 
 
    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)
    
    # push!(e, ei)
    # push!(angles, α)
    println(length(generators))
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei), " Error: ",real(m[1]) -  real(ei) )
    return real(ei)
end

molecule()
ansatz()
