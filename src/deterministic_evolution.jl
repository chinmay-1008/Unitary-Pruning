using Random, Distributions
using BenchmarkTools
using LinearAlgebra
using PauliOperators
"""
    deterministic_pauli_rotations(generators::Vector{P}, angles, o::P, ket; nsamples=1000) where {N, P<:Pauli}


"""

@inline commute(p1::Union{Pauli{N}, PauliBasis{N}}, p2::Union{Pauli{N}, PauliBasis{N}}) where {N} = iseven(count_ones(p1.x & p2.z) - count_ones(p1.z & p2.x)) 

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


function newclip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(p->(abs(p.second) ≥ thresh) , ps)
end

function bfs_evolution(generators::Union{Vector{Pauli{N}}, Vector{PauliBasis{N}}}, angles, o::PauliSum{N}, ket ; thresh=1e-3) where {N}

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

        for (oi,coeff) in o_transformed
           
            abs(coeff) > thresh || continue

            if commute(oi, PauliBasis(g)) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
            temp_norm2+=abs(coeff)^2
        end
        sum!(o_transformed, sin_branch) 
        newclip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)
    # println(coeff_norm2)
    # println("dfhvskjdvjksdfsdjkfsjhfkshfleshfkehfkl")

    return expval, n_ops, coeff_norm2
end

function weight(ps::Union{Pauli{N}, PauliBasis{N}}) where {N}
    x = ps.x
    z = ps.z
    return count_ones(x | z)
end

function majorana_weight(Pb::Union{Pauli{N}, PauliBasis{N}}) where N

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


function weightclip!(ps::PauliSum{N}; lc = 0) where {N}
    filter!(p-> weight(p.first) ≤ lc , ps)
end

function majorana_clip!(ps::PauliSum{N}; lc = 0) where {N}
    filter!(p-> majorana_weight(p.first) ≤ lc , ps)
end

function myclip!(ps::PauliSum{N}; thresh=1e-16, lc = 0, w_type = 0) where {N}
    if w_type == 0 
        filter!(p->(weight(p.first) ≤ lc) && (abs(p.second) ≥ thresh) , ps)
    else
        filter!(p->(majorana_weight(p.first) ≤ lc) && (abs(p.second) ≥ thresh) , ps)
    end     
end
# w_type is the type of clipping, 0 is Pauli weight clipping and 1 is majorana weight clipping
function bfs_evolution_weight(generators::Union{Vector{Pauli{N}}, Vector{PauliBasis{N}}}, angles, o::PauliSum{N}, ket ; thresh=1e-3, w_type = 0, w = 2) where {N}

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

        for (oi,coeff) in o_transformed
           
            abs(coeff) > thresh || continue

            if commute(oi, g) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos[t]

                # sin branch
                oj = g * oi    # multiply the pauli's
                sum!(sin_branch, oj * vsin[t] * coeff * 1im)
  
            end
        end
        sum!(o_transformed, sin_branch) 
        # clip!(o_transformed, thresh=thresh)
        # myclip!(o_transformed, thresh=thresh, lc = w, w_type = w_type)
        if w_type == 0

            weightclip!(o_transformed, lc = w)

        elseif w_type == 1 

            majorana_clip!(o_transformed, lc = w)

        end 

        n_ops[t] = length(o_transformed)
    end

    coeff_norm2 = 0

    for (oi,coeff) in o_transformed
        expval += coeff*expectation_value(oi, ket)
        coeff_norm2+= abs(coeff)^2      # final list of operators
    end
    coeff_norm2 = sqrt(coeff_norm2)

    return expval, n_ops, coeff_norm2
end
