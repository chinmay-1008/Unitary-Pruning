using Random, Distributions
using BenchmarkTools
using LinearAlgebra
"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}


"""
function deterministic_pauli_rotations(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thres=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    
    vcos = cos.(angles)
    vsin = sin.(angles)

    # collect our results here...
    expval = zero(ComplexF64)


    # opers = PauliSum{N}(Dict(o.pauli=>1.0*(1im)^o.θ))#Dict{Tuple{Int128, Int128}, Complex{Float64}}((o.pauli.z,o.pauli.x)=>1.0*(1im)^o.θ)
    opers = PauliSum(o)
  
    n_ops = zeros(Int,nt)
    
    for t in 1:nt

        g = generators[t]
        branch_opers = PauliSum(N)

        sizehint!(branch_opers.ops, 1000)
        for (key,value) in opers.ops
            
            oi = key
            if commute(oi, g.pauli)
                if haskey(branch_opers, oi)
                    branch_opers[oi] += value
                else
                    branch_opers[oi] = value
                end
                continue
            end
            if abs(value) > thres #If greater than threshold then split the branches
                # cos branch
                coeff = vcos[t] * value
                           
                if haskey(branch_opers, oi) # Add operator to dictionary if the key doesn't exist
                    branch_opers[oi] += coeff
                else
                    branch_opers[oi] = coeff #Modify the coeff if the key exists already
                end

                # sin branch
                coeff = vsin[t] * 1im * value
                oi = Pauli{N}(0, oi)
                oi = g * oi    # multiply the pauli's

                if haskey(branch_opers, oi.pauli) # Add operator to dictionary if the key doesn't exist
                    branch_opers[oi.pauli] += coeff * (1im)^oi.θ
                else
                    branch_opers[oi.pauli] = coeff * (1im)^oi.θ #Modify the coeff if the key exists already
                end
            end
        end
        n_ops[t] = length(branch_opers)
        opers = deepcopy(branch_opers) # Change the list of operators to the next row

    end

    for (key,value) in opers.ops
        oper = Pauli(UInt8(0), key)
        expval += value*PauliOperators.expectation_value(oper, ket)
    end
   
    return expval, n_ops
end




"""
    bfs_evolution(generators::Vector{Pauli{N}}, angles, o::Pauli{N}, ket ; thres=1e-3) where {N}


"""
function bfs_evolution(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3) where {N}

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
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2
end

# here layer = false means the threshold is set randomly for every pauli in every layer, when true threshold is randomly set for every layer
function stochastic_bfs(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket; sigma=1e-3, layer::Bool = false) where {N}

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
        

            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end

        sum!(o_transformed, sin_branch)

        if layer
            # similar to threshold but based on a distribution
            rn = abs(rand(Normal(0, sigma)))
            clip!(o_transformed, thresh=rn)
        else
    
            # killing random paulis based on the coefficients 
            for (oi, coeff) in o_transformed.ops
                rn = abs(rand(Normal(0, sigma)))

                if abs(coeff) < rn

                    delete!(o_transformed.ops, oi)

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


# doesn't work as intended
function flat_prob_kill_bfs(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket; threshold=1e-3) where {N}

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
        

            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end

        sum!(o_transformed, sin_branch)

            # killing random paulis based on the coefficients 
        for (oi, coeff) in o_transformed.ops
            rn = rand()

            if rn < threshold

                delete!(o_transformed.ops, oi)

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


function stochastic_bfs_flat(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket; sigma=1e-3, layer::Bool = false) where {N}

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
        

            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end

        sum!(o_transformed, sin_branch)

        if layer
            # similar to threshold but based on a distribution
            rn = abs(rand(Uniform(0, sigma)))
            clip!(o_transformed, thresh=rn)
        else
    
            # killing random paulis based on the coefficients 
            for (oi, coeff) in o_transformed.ops
                rn = abs(rand(Uniform(0, sigma)))

                if abs(coeff) < rn

                    delete!(o_transformed.ops, oi)

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


function stochastic_bfs_renormalized(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket; sigma=1e-3, layer::Bool = false) where {N}

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
        

            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end

        sum!(o_transformed, sin_branch)
        temp_coeff_norm2 = 0

        if layer
            # similar to threshold but based on a distribution
            rn = abs(rand(Uniform(0, sigma)))
            clip!(o_transformed, thresh=rn)
        else
    
            # killing random paulis based on the coefficients 
            for (oi, coeff) in o_transformed.ops
                rn = abs(rand(Uniform(0, sigma)))

                if abs(coeff) < rn

                    delete!(o_transformed.ops, oi)

                else

                    temp_coeff_norm2+= abs(coeff)^2   

                end
            end
        end
        temp_coeff_norm2 = sqrt(temp_coeff_norm2)

        for (oi, coeff) in o_transformed.ops
            o_transformed[oi] = coeff/temp_coeff_norm2
        
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



function stochastic_bfs_coeff_error(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket; sigma=1e-3, epsilon= 1e-3, layer::Bool = false) where {N}

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
        

            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end

        sum!(o_transformed, sin_branch)

        for (oi , coeff) in o_transformed.ops
            rn_c = rand(Normal(0, epsilon))
            temp_coeff = 0
            if real(coeff) == 0
                temp_coeff += 1im*rn_c
            elseif imag(coeff) == 0
                temp_coeff += rn_c
            else
                temp_coeffcoeff += rnc + 1im*rn_c
            end

            o_transformed[oi] = coeff + temp_coeff
        end

        if layer
            # similar to threshold but based on a distribution
            rn = abs(rand(Normal(0, sigma)))
            clip!(o_transformed, thresh=rn)
        else
    
            # killing random paulis based on the coefficients 
            for (oi, coeff) in o_transformed.ops

                rn = abs(rand(Normal(0, sigma)))


                if abs(coeff) <  rn

                    delete!(o_transformed.ops, oi)

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


function bfs_coeff_error(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3, epsilon = 1e-3) where {N}

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
        temp_norm2 = 0

        for (oi,coeff) in o_transformed.ops
           
            abs(coeff) > thresh || continue


            if commute(oi, g.pauli) == false

                rn_c = rand(Normal(0, epsilon))

                # if real(coeff) == 0
                #     coeff += 1im*rn_c
                # elseif imag(coeff) == 0
                #     coeff += rn_c
                # else
                #     coeff += rnc + 1im*rn_c
                # end

                # cos branch
                o_transformed[oi] = coeff * (vcos[t] + rn_c)

                # sin branch
                oj = g * oi    # multiply the pauli's
                rn_s = rand(Normal(0, epsilon))

                sum!(sin_branch, oj * (vsin[t]+ rn_s) * coeff * 1im )
  
            end
            temp_norm2+=abs(coeff)^2
        end


        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2
end


function bfs_coeff_error_renormalized(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3, epsilon = 1e-3) where {N}

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
        temp_norm2 = 0

        for (oi,coeff) in o_transformed.ops
           
            abs(coeff) > thresh || continue


            if commute(oi, g.pauli) == false

                rn_c = rand(Normal(0, epsilon))

                # if real(coeff) == 0
                #     coeff += 1im*rn_c
                # elseif imag(coeff) == 0
                #     coeff += rn_c
                # else
                #     coeff += rnc + 1im*rn_c
                # end

                # cos branch
                o_transformed[oi] = coeff * (vcos[t] + rn_c)

                # sin branch
                oj = g * oi    # multiply the pauli's
                rn_s = rand(Normal(0, epsilon))

                sum!(sin_branch, oj * (vsin[t]+ rn_s) * coeff * 1im )
  
            end
        end


        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)

        for (oi, coeff) in o_transformed.ops
            temp_norm2 += abs(coeff)^2
        end

        temp_norm2 = sqrt(temp_norm2)

        for (oi , coeff) in o_transformed.ops
            o_transformed[oi] = coeff/temp_norm2
        end
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2
end



function bfs_angle_error(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3, epsilon = 1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)

    for i in 1:nt
        rn_a = rand(Normal(0, epsilon))
        angles[i] += rn_a
    end
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
        temp_norm2 = 0

        for (oi,coeff) in o_transformed.ops
           
            abs(coeff) > thresh || continue


            if commute(oi, g.pauli) == false

                rn_c = rand(Normal(0, epsilon))

                # if real(coeff) == 0
                #     coeff += 1im*rn_c
                # elseif imag(coeff) == 0
                #     coeff += rn_c
                # else
                #     coeff += rnc + 1im*rn_c
                # end

                # cos branch
                o_transformed[oi] = coeff * (vcos[t])

                # sin branch
                oj = g * oi    # multiply the pauli's
                rn_s = rand(Normal(0, epsilon))

                sum!(sin_branch, oj * (vsin[t]) * coeff * 1im )
  
            end
            temp_norm2+=abs(coeff)^2
        end


        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2
end

function weight(ps::FixedPhasePauli{N}) where {N}
    x = ps.x
    z = ps.z
    return count_ones(x | z)
end

function myclip!(ps::PauliSum{N}; thresh=1e-16, lc = 0) where {N}
    filter!(p->(abs(p.second) > thresh) | (weight(p.first) < lc), ps.ops)
end


function bfs_evolution_diff(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3, γ = 0, lc = 0) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    vcos = cos.(angles)
    vsin = sin.(angles)
    # a = rand(3)*5
    # collect our results here...
    expval = zero(ComplexF64)


    o_transformed = deepcopy(o)
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
        
        for (oi, coeff) in o_transformed.ops
            # x = oi.x
            # z = oi.z
            # ls = count_ones(x | z)
            ls = weight(oi)
            if ls > lc
                # coeff *= exp(-γ * (ls - lc))
                # if abs(coeff) < thresh
                #     delete!(o_transformed.ops, oi)
                # else
                o_transformed[oi] = coeff * exp(-γ * (ls - lc))
                # end
            end
        end
        # clip!(o_transformed, thresh=thresh)
        myclip!(o_transformed, thresh = thresh, lc=lc)
        n_ops[t] = length(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2
end

function dissipation(o::PauliSum{N}; γ = 0.1, lc = 0) where N
    o_transformed_d = deepcopy(o)
    # temp_coeff = Float64(0)
    for (oi_d, coeff_d) in o_transformed_d.ops
        ls = weight(oi_d)
        if ls >= lc
            # temp_coeff = coeff_d * exp(-γ*(ls - lc))
            o_transformed_d[oi_d] = coeff_d * exp(-γ*(ls - lc))
            # println(exp(-γ*(ls - lc)))
            # println("Transform: ", coeff_d, " to: ", coeff_d * exp(-γ*(ls - lc)))
        end
    end 

    return o_transformed_d
end


function bfs_evolution_new_diff(generators::Vector{Pauli{N}}, angles, o, hj, k; thresh=1e-3,γ = 0, lc = 0) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    vcos = cos.(2.0*angles)
    vsin = sin.(2.0*angles)
 
    initial_state = deepcopy(o)
    sin_branch = PauliSum(N)
 
    n_ops = zeros(Int,nt)
    one_step = 4*N
    cj = []
    t_cj = 0.0
    temp_norm = 0

    for ki in 1:k
        println("Step: ", ki, " of: ", k)

        if ki == 1
            for h in hj
                state = adjoint(h) * initial_state
                c_i = tr(state)
                temp_norm += abs(c_i)
            end
        end

        temp_cj = []
        o_step = PauliSum(N)

        for (oi_state, coeff_state) in initial_state.ops
            o_transformed = PauliSum(N)
            o_transformed += oi_state
            # display(o_transformed)

            for step in 1:one_step

                g = generators[step]
    
                sin_branch = PauliSum(N)
    
                for (oi,coeff) in o_transformed.ops
                    # abs(coeff) > thresh || continue
    
    
                    if commute(oi, g.pauli) == false
                        
                        # cos branch
                        o_transformed[oi] = coeff * vcos[step]
    
                        # sin branch
                        oj = g * oi    # multiply the pauli's
                        sum!(sin_branch, oj * coeff * vsin[step] * 1im)
        
                    end
                end
                sum!(o_transformed, sin_branch)
                
                # clip!(o_transformed, thresh=thresh)
                # myclip!(o_transformed, thresh = thresh, lc=lc)
            end
            o_transformed = dissipation(o_transformed, γ = γ, lc = lc)
            o_transformed = mul!(o_transformed, coeff_state)
            o_step += o_transformed

        end 
    
        # println("Calculating the trace for C_j")
        for h in hj
            t_cj = tr(adjoint(h)*o_step)
            push!(temp_cj, t_cj/temp_norm)
        end
          
        initial_state = deepcopy(o_step)
        println(temp_cj)
        push!(cj, temp_cj)

        println("--------------------------------------------------")

    end

    return cj
end


    # for t in 1:nt

    #     g = generators[t]

    #     sin_branch = PauliSum(N)

    #     for (oi,coeff) in o_transformed.ops
           
    #         abs(coeff) > thresh || continue


    #         if commute(oi, g.pauli) == false
                
    #             # cos branch
    #             o_transformed[oi] = coeff_i * coeff * vcos[t]

    #             # sin branch
    #             oj = g * oi    # multiply the pauli's
    #             sum!(sin_branch, oj * coeff_i * vsin[t] * coeff * 1im)
  
    #         end
    #     end
    #     sum!(o_transformed, sin_branch)

    #     if t%(4*N) == 0
    #         # println("Dissipating")
    #         # println(t)
    #         o_transformed = dissipation(o_transformed, γ = γ, lc = lc)
    #     end
        
    #     # clip!(o_transformed, thresh=thresh)
    #     # myclip!(o_transformed, thresh = thresh, lc=lc)
    #     n_ops[t] = length(o_transformed)
    # end
    # Final timestep without dissipation
    # for t in 1:nt

    #     g = generators[t]

    #     sin_branch = PauliSum(N)

    #     for (oi,coeff) in o_transformed.ops
           
    #         abs(coeff) > thresh || continue


    #         if commute(oi, g.pauli) == false
                
    #             # cos branch
    #             o_transformed[oi] = coeff_i * coeff * vcos[t]

    #             # sin branch
    #             oj = g * oi    # multiply the pauli's
    #             sum!(sin_branch, oj * coeff_i * vsin[t] * coeff * 1im)
  
    #         end
    #     end
    #     sum!(o_transformed, sin_branch)

    #     if t%(4*N) == 0
    #         # o_transformed = dissipation(o_transformed, γ = γ, lc = lc)
    #         n_ops[t] = length(o_transformed) 
    #         break
    #     end
        
    #     # clip!(o_transformed, thresh=thresh)
    #     # myclip!(o_transformed, thresh = thresh, lc=lc)
    #     n_ops[t] = length(o_transformed)
    # end

    # println(coeff_norm2)


