function bfs_subspace_expansion(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3) where {N}
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
        n_ops[t] = length(o_transformed)
    end
    temp_norm = 0
    for (oi, coeff) in o_transformed.ops
        temp_norm += abs(coeff)^2
    end
    temp_norm = sqrt(temp_norm)

#adding a subspace evolution
    for t in 1:nt
        g = generators[t]

        for (oi,coeff) in o_transformed.ops
            coeff = coeff/temp_norm
            if commute(oi, g.pauli)
                continue
            else
                temp = g * oi
                o_transformed[oi] = (o_transformed[oi]/temp_norm) + coeff * vcos[t]
                # println(o_transformed[oi])
                if haskey(o_transformed.ops, temp.pauli)    
                    o_transformed[temp.pauli] = (get(o_transformed, temp.pauli)/temp_norm) +  1im*vsin[t]*coeff*get_phase(temp)*2
                    println(o_transformed[temp.pauli])
                # else
                    # println("sin")
                    # o_transformed[temp] = 1im * vsin[t] * coeff

                end
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