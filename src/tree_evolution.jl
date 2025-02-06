using Random, Distributions
using BenchmarkTools
using LinearAlgebra

function tree_evolution(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3) where {N}

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
    dicts = Vector(undef, nt)    
    sin_branch = PauliSum(N)
 
    n_ops = zeros(Int,nt)
    
    for t in 1:nt

        g = generators[t]

        sin_branch = PauliSum(N)
        temp_norm2 = 0

        for (oi,coeff) in o_transformed.ops
           
            abs(coeff) > thresh || continue


            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
            temp_norm2+=abs(coeff)^2
        end
        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
        println((n_ops[t]))
        # println(o_transformed)
        dicts[t] = deepcopy(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2, dicts
end


function tree_evolution_otoc(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    vcos = cos.(angles)
    vsin = sin.(angles)

    # collect our results here...
    step_evo = Int64(nt/(3*(N-1)))
    otoc = Vector{Float64}(undef, step_evo)
    temp_o = PauliSum(N)
    time = 1
    o_transformed = deepcopy(o)
    dicts = Vector(undef, nt)    
    sin_branch = PauliSum(N)
 
    
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
        # println((n_ops[t]))
        # println(o_transformed)
        dicts[t] = deepcopy(o_transformed)
        println("Generator: ", t)
        if t%(3*(N-1)) == 0
            println(time)
            os = o*o_transformed
            println("First otoc done")
            # display(os)
            temp_o = os*os
            println("Second otoc done")
            otoc[time] = real(expectation_value(temp_o, ket))
            time += 1
        end

    end

    return otoc, dicts
end