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

function otoc_expectation_value(o_oto_1::PauliSum{N}, o_otoc_2::PauliSum{N}, ket) where {N}
    output = 0

    for (oi_1, coeff_1) in o_oto_1.ops
        for (oi_2, coeff_2) in o_otoc_2.ops
            output += expectation_value(oi_1*oi_2, ket)*coeff_1*coeff_2
        end
    end

    return output
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
    otoc = Vector{Float64}(undef, step_evo+1)
    temp_o = PauliSum(N)
    time = 1+1
    o_transformed = deepcopy(o)
    dicts = Vector(undef, nt)    
    sin_branch = PauliSum(N)
    
    otoc[1] = real(expectation_value(o, ket))
    
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
        dicts[t] = deepcopy(o_transformed)

        if t%(3*(N-1)) == 0
            println("Time: ", time)
            os = o*o_transformed
            println("First otoc done")
            # display(os)
            # temp_o = os*os
            println("Second otoc done")
            # otoc[time] = real(expectation_value(temp_o, ket))
            otoc[time] = real(otoc_expectation_value(os, os, ket))
            # otoc[time] = tr(adjoint(temp_o)*temp_o)

            time += 1
        end

    end

    return otoc, dicts
end