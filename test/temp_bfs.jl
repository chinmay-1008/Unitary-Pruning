using UnitaryPruning
using Random, Distributions
using PauliOperators
using Printf
using Plots


function bfs_new(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket , thresh , sigma) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    vcos = cos.(angles)
    vsin = sin.(angles)

    # collect our results here...
    expval = zero(ComplexF64)


    o_transformed = deepcopy(o)
    sin_branch = PauliSum(N)
 
    n_ops = zeros(Int,nt)
    
    for t in 1:nt

        g = generators[t]

        sin_branch = PauliSum(N)

        for (oi,coeff) in o_transformed.ops
           
            abs(coeff) > thresh || continue


            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end

        sum!(o_transformed, sin_branch)
        clip!(o_transformed, thresh=thresh)

        for (oi, coeff) in o_transformed.ops
            rn = abs(rand(Normal(0, sigma)))

            if abs(coeff) < rn
                o_transformed[oi] = 0

            end
        end


        n_ops[t] = length(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)

    return expval, n_ops, coeff_norm2
end

function run(; N=10, k=5, thresh=1e-3, sigma = 1e-3)
   
    ket = KetBitString(N, 0) 
    o = Pauli(N, Z=[1])

    angles = [] 
    e = [] 
    
    i = 4
    α = i * π / 32 
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = α, Jy = α,Jz = α, k=k)
    # ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=thresh)

    ei_n , nops_n, c_norm2_n = bfs_new(generators, parameters, PauliSum(o), ket, thresh, sigma)
    
    # push!(e, ei)
    # push!(angles, α)
    # @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.5f\n", α, real(ei), imag(ei), maximum(nops), c_norm2, thresh)
    @printf(" α: %6.4f e: %12.8f+%12.8fi nops: %6i norm2: %3.8f threshold: %3.10f\n", α, real(ei_n), imag(ei_n), maximum(nops_n), c_norm2_n, thresh)

    # U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    # o_mat = Matrix(o)
    # m = diag(U'*o_mat*U)
    # println(real(m[1]), " ", real(ei), " Error: ",real(m[1]) -  real(ei) )
    return real(ei_n)
end

function run_temp()
    stds = [i*1e-3 for i in 1:20]
    energy = []
    for sigma in stds

        ei = run(N = 10, k = 10, thresh = 1e-3, sigma=sigma)
        push!(energy, ei)
    end

    println(energy)
    scatter(stds, energy, smooth=true)
    xlabel!("σ")
    ylabel!("Expectation Value")
    title!("E{Z_1} vs σ (N=10; k=10; α=4π/32)")
    savefig("test/ener_vs_sigma_hei.pdf")

end

run_temp()