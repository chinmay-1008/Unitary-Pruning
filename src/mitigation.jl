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
    
    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
    end

    expval_bfs = expval
    # display(o_transformed)
    # @printf(" Expectation value for BFS: %12.8f\n", expval)

    expval = 0
    # println(n_ops[end])

    # Now go back to t=0
    for (oi,coeff) in o_transformed.ops
        o_transformed[oi] = 0
    end
   
    sum!(o_transformed, o)
    
    #
    #adding a subspace evolution
    for t in 1:nt
        g = generators[t]

        next_layer = PauliSum(N)
        for (oi,coeff) in o_transformed.ops
            if commute(oi, g.pauli) == false
                temp = g * oi
                o_transformed[oi] = coeff * vcos[t]
                if haskey(o_transformed.ops, temp.pauli)
                    # o_transformed[temp.pauli] = get(o_transformed, temp.pauli) +  1im*vsin[t]*coeff*get_phase(temp)
                    sum!(next_layer, temp * vsin[t] * coeff * 1im)
                end
            end
        end
        sum!(o_transformed, next_layer)

        n_ops[t] = length(o_transformed)
    end
    # println(n_ops[end])

    # display(o_transformed)
    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)


    return expval, expval_bfs, n_ops, coeff_norm2


end