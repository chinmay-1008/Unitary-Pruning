using Random, Distributions
using BenchmarkTools
using LinearAlgebra
using PauliOperators
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

function weight(ps::FixedPhasePauli{N}) where {N}
    x = ps.x
    z = ps.z
    return count_ones(x | z)
end

function majorana_weight(Pb::Union{Pauli{N}, FixedPhasePauli{N}}) where N

    w = 0
    control = true
    # tmp = Pb.z & ~Pb.x  # Bitwise AND with bitwise NOT
    Ibits = ~(Pb.z|Pb.x)
    Zbits = Pb.z & ~Pb.x
    for i in reverse(1:N)  # Iterate from N down to 1
        xbit = (Pb.x >> (i - 1)) & 1 != 0
        Zbit = (Zbits >> (i - 1)) & 1 != 0
        Ibit = (Ibits >> (i - 1)) & 1 != 0
        #println("i=$i, xbit=$xbit, Zbit=$Zbit, Ibit=$Ibit, control=$control, w=$w")
        if Zbit && control || Ibit && !control
            w += 2
        elseif xbit
            control = !control
            w += 1
        end
    end
    return w
end

function myclip!(ps::PauliSum{N}; thresh=1e-16, lc = 0, w_type = 0) where {N}
    if w_type == 0 
        filter!(p->(weight(p.first) ≤ lc) && (abs(p.second) ≥ thresh) , ps.ops)
    else
        filter!(p->(majorana_weight(p.first) ≤ lc) && (abs(p.second) ≥ thresh) , ps.ops)
    end     
end

function weightclip!(ps::PauliSum{N}; lc = 0) where {N}
    filter!(p-> weight(p.first) <= lc , ps.ops)
end

function majorana_clip!(ps::PauliSum{N}; lc = 0) where {N}
    filter!(p-> majorana_weight(p.first) <= lc , ps.ops)
end

# w_type is the type of clipping, 0 is Pauli weight clipping and 1 is majorana weight clipping
function bfs_evolution_weight(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket ; thresh=1e-3, w_type = 0, w = 2) where {N}

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
        # println("Initial ", length(o_transformed))
        # display(o_transformed)
        # clip!(o_transformed, thresh=thresh)
        myclip!(o_transformed, thresh=thresh, lc = w, w_type = w_type)
        # if w_type == 0

        #     weightclip!(o_transformed, lc = w)

        # elseif w_type == 1 

        #     majorana_clip!(o_transformed, lc = w)

        # end

        # println("Clipped ", length(o_transformed))
        # display(o_transformed)
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
    for (oi_d, coeff_d) in o.ops
        ls = weight(oi_d)
        if ls >= lc
            o[oi_d] = coeff_d* exp(-γ*(ls - lc))
        end
    end 
    return o
end


function bfs_evolution_new_diff(generators::Vector{Pauli{N}}, angles, initial_state, hj, k, dt, T; thresh=1e-3,γ = 0, lc = 0) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    vcos = cos.(2.0*angles)
    vsin = sin.(2.0*angles)
  
    one_step = 4*N
    cj = Vector{Vector{Float64}}(undef, Int64(k))
    temp_norm = 0

    if !isempty(hj)
        for h in hj
            state = adjoint(h) * initial_state
            temp_norm += abs(tr(state))
        end
    end

    for ki in 1:k
        println("Step: ", ki, " of: ", k)
        # println("Number of ops: ", length(initial_state))

        temp_cj = Vector{Float64}(undef, length(hj))

        o_step = PauliSum(N)
        l = ReentrantLock()
        
        # Threads.@threads for (oi_state, coeff_state) in collect(initial_state.ops)
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
                
                # myclip!(o_transformed, thresh = thresh, lc=lc)
            end
            mul!(o_transformed, coeff_state)
            o_transformed = dissipation(o_transformed, γ = γ, lc = lc)

            clip!(o_transformed, thresh=thresh)

            # lock(l) do
            #     sum!(o_step, o_transformed)
            # end
            sum!(o_step, o_transformed)

        end 
    
        for (idx, h) in enumerate(hj)
            temp_cj[idx] = tr(adjoint(h) * o_step) / temp_norm
        end
          
        initial_state = o_step
        cj[Int64(ki)] = temp_cj

        println("--------------------------------------------------")

    end
    return cj, initial_state
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

function dot_pdt(bra_state, ket_state)
    val = zero(ComplexF64)
    for (ket, coeff_ket) in ket_state
        for (bra, coeff_bra) in bra_state
            if ket.v == bra.v
                val += (coeff_bra' * coeff_ket)
            end
        end
    end
    return val
end 

function bfs_evolution_central(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket, dt, pi_pulse ; thresh=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    angles *= (2*dt)
    
    for j in pi_pulse
        # println(angles[j])
        angles[j] /= dt
        # println(generators[j], " ", angles[j])
        # println(sin(angles[j]))
    end

    vcos = cos.(angles)
    vsin = sin.(angles)


    # collect our results here...
    expval = zero(ComplexF64)

    o_transformed = deepcopy(o)
    sin_branch = PauliSum(N)
    d_nt = Int64(floor((nt-length(pi_pulse))/(3*(N-1))))
    # println(d_nt)
    dicts = Vector(undef, d_nt)    
    n_ops = zeros(Int,nt)
    ham_ops = 0
    exp_time = Vector{Float64}([])

    temp_state = o_transformed * ket
    exp_l = dot_pdt(ket, temp_state)
    display(temp_state)
    println("Expectation Valuel: ", exp_l)
    push!(exp_time, abs(exp_l))

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
                oj = oi * g    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
    
            end
            temp_norm2+=abs(coeff)^2
        end
        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
        ham_ops +=1

        if ham_ops % (3*(N-1)) == 0
            # println(ham_ops)
            # println("op: ", g, " t: ", t, " ham: ", ham_ops)
            temp_state = o_transformed * ket

            exp_l = dot_pdt(ket, temp_state)
            # exp_l = 0
            # for (state, s_coeff) in ket
                # exp_l += (abs(s_coeff)^2) * expectation_value(o_transformed, state) 
            # end
            # exp_l = expectation_value(o_transformed, ket)
            # println("--------------------taking exp-------------------")
            push!(exp_time, abs(exp_l))
            # display(temp_state)
            # println(abs(exp_l))
            dicts[Int(ham_ops/(3*(N-1)))] = deepcopy(o_transformed)
        end

        if g == Pauli(N, X = [1])
            # println("-----------------------------")
            println(g)
            ham_ops -= 1
            println(ham_ops)
            println("-----------------------------")
            # exp_l = 0
            # for (state, s_coeff) in ket
                # exp_l += (abs(s_coeff)^2) * expectation_value(o_transformed, state) 
            # end
            # println("Pulse Exp Value: ", exp_l)
        end
        # println(t, " ", ham_ops)
        # println("opes done: ", g)
        # display(o_transformed)
        # println("-----------------------------")

    end


    # for (state, s_coeff) in ket
    #     display(state)
    #     coeff_norm2 = 0
    #     temp_exp = 0
    #     for (oi,coeff) in o_transformed.ops
    #         temp_exp += coeff*expectation_value(oi, state)
    #         coeff_norm2+= abs(coeff)^2      # final list of operators
    #     end 
    #     temp_exp *= s_coeff^2
    #     expval += temp_exp
    # end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
    #     expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end 
    temp_state = o_transformed * ket
    # exp_state = o * temp_state
    # display(temp_state)
    # println("Exp: ", dot_pdt(ket, temp_state))
    # println("Purity S+: ", abs(dot_pdt(temp_state, temp_state)))

    expval = dot_pdt(ket, temp_state)
    coeff_norm2 = sqrt(coeff_norm2)

    return expval, n_ops, coeff_norm2, exp_time, dicts
end

function bfs_evolution_central_state(generators::Vector{Pauli{N}}, angles, o::PauliSum{N}, ket, dt, pi_pulse ; thresh=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
    # reverse!(generators)
    # reverse!(angles)
    nt = length(angles)
    length(angles) == nt || throw(DimensionMismatch)
    angles *= (dt)
    
    for j in pi_pulse
        # println(angles[j])
        angles[j] /= dt
        # println(generators[j], " ", angles[j])
        # println(sin(angles[j]))
    end

    vcos = cos.(angles)
    vsin = sin.(angles)


    # collect our results here...
    expval = zero(ComplexF64)

    o_transformed = PauliSum(N) + Pauli(N)
 
    sin_branch = PauliSum(N)
    d_nt = Int64(floor((nt-length(pi_pulse))/(3*(N-1))))
    # println(d_nt)
    dicts = Vector(undef, d_nt)    
    n_ops = zeros(Int,nt)
    ham_ops = 0
    exp_time = Vector{Float64}([])

    temp_state = o_transformed * ket

    exp_state = o * temp_state
    exp_l = dot_pdt(temp_state, exp_state)
    display(exp_state)
    println("Exp Val: ", abs(exp_l))
    push!(exp_time, abs(exp_l))

    for t in 1:nt
        
        g = generators[t]

        sin_branch = PauliSum(N)
        temp_norm2 = 0

        for (oi,coeff) in o_transformed.ops
            
            abs(coeff) > thresh || continue


            # if commute(oi, g.pauli) == false
                
                # cos branch
            o_transformed[oi] = coeff * vcos[t]

                # sin branch
            oj = oi * g    # multiply the pauli's
            sum!(sin_branch, oj * vsin[t] * coeff * (-1im))
    
            # end
            temp_norm2+=abs(coeff)^2
        end
        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
        ham_ops +=1

        if ham_ops % (3*(N-1)) == 0
            # println(ham_ops)
            # println("op: ", g, " t: ", t, " ham: ", ham_ops)
            temp_state = o_transformed * ket
            # display(temp_state)
            # println("--------------------taking exp-------------------")

            exp_state = o * temp_state
            exp_l = dot_pdt(temp_state, exp_state)

            push!(exp_time, abs(exp_l))
            # println(abs(exp_l))
            dicts[Int(ham_ops/(3*(N-1)))] = deepcopy(o_transformed)
        end

        if g == Pauli(N, X = [1])
            println("-----------------------------")
            println(g)
            ham_ops -= 1
            println(t, " ", ham_ops)
            println("-----------------------------")
        end

    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed.ops
        # expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end 

    temp_state = o_transformed * ket
    exp_state = o * temp_state
    display(temp_state)
    println("Purity: ", abs(dot_pdt(temp_state, temp_state))^2)
    println("Purity S+: ", abs(dot_pdt(exp_state, exp_state)))

    expval = dot_pdt(temp_state, exp_state)
    coeff_norm2 = sqrt(coeff_norm2)

    return expval, n_ops, coeff_norm2, exp_time, dicts
end