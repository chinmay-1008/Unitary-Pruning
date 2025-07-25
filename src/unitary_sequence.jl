function get_unitary_sequence_1D(o::Pauli{N}; α=.01, k=10) where N
    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])
    # print("alpha", α, "\n")
    # Loop over trotter steps
    for ki in 1:k
        ## ZZ layer
        # e^{i π/2 P2} e^{i π P1 /2}|ψ>
        for i in 1:N-1
            pi = Pauli(N, Z=[i, i + 1])
            push!(generators, pi)
            push!(parameters, π/2)
        end
        #pbc 
        pi = Pauli(N, Z=[N, 1])
        push!(generators, pi)
        push!(parameters, π/2)

        ## X layer
        # e^{i αn (-X) / 2}
        for i in 1:N
            pi = Pauli(N, X=[i])
            pi = Pauli{N}((pi.θ + 2)%4, pi.pauli) # this accounts for the fact that the papers have -X and positive ZZ
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    return generators, parameters
end

function kicked_ising_2D(L::Int; α = 0.01, k = 10)
    N = L * L  # Total number of spins
    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])

    # Map 2D coordinates to 1D index
    index(i, j) = ((i - 1) % L) * L + ((j - 1) % L) + 1

    for step in 1:k
        ## ZZ layer 
        for i in 1:L, j in 1:L
            # Horizontal (right neighbor)
            ni = index(i, j)
            nj = index(i, j + 1)
            push!(generators, Pauli(N, Z=[ni, nj]))
            push!(parameters, π / 2)

            # Vertical (bottom neighbor)
            ni = index(i, j)
            nj = index(i + 1, j)
            push!(generators, Pauli(N, Z=[ni, nj]))
            push!(parameters, π / 2)
        end

        ## X layer 
        for i in 1:N
            pi = Pauli(N, X=[i])
            pi = Pauli{N}((pi.θ + 2) % 4, pi.pauli)  # handles -X
            push!(generators, pi)
            push!(parameters, α)
        end
    end

    return generators, parameters
end



function get_unitary_sequence_2D(o::Pauli{N}; α=.01, k=10) where N


    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])
    bridges = [[1,5,7], [3,6,9]]
    sequences = [[1,3], [7,9]]

    for ki in 1:k
        # e^{i π/4 P2} e^{i π P1 /2}|ψ>
        ## ZZ layer
        for qubit in sequences
            for i in qubit[1]:qubit[2]
                pi = Pauli(N, Z=[i, i + 1])
                push!(generators, pi)
                push!(parameters, π/2)
            end
        end
        #bridges
        for link in bridges
            pi = Pauli(N, Z=[link[1], link[2]])
            push!(generators, pi)
            push!(parameters, π/2)
            
            pi = Pauli(N, Z=[link[2], link[3]])
            push!(generators, pi)
            push!(parameters, π/2)
        end
        ## X layer
        # e^{i αn Pn / 2}
        for i in 1:N
            pi = Pauli(N, X=[i])
            pi = Pauli{N}((pi.θ + 2)%4, pi.pauli) # this accounts for the fact that the papers have -X and positive ZZ
            push!(generators, pi)
            push!(parameters, α)
        end
    end
    return generators, parameters
end



function build_time_evolution_matrix(generators::Vector{Pauli{N}}, angles::Vector) where N
    U = Matrix(Pauli(N))
    nt = length(generators)
    length(angles) == nt || throw(DimensionMismatch)
    for t in 1:nt
        α = angles[t]
        U = cos(α/2) .* U .- 1im*sin(α/2) .* U * Matrix(generators[t])

    end

    return U 
end

function build_time_evolution_matrix_fast(generators::Vector{Pauli{N}}, angles::Vector{<:Real}) where N
    nt = length(generators)
    length(angles) == nt || throw(DimensionMismatch())

    # Start with identity matrix of size 2^N × 2^N
    U = Matrix{ComplexF64}(I, 2^N, 2^N)
    W = Matrix{ComplexF64}(undef, 2^N, 2^N)  # workspace to avoid allocs

    for t in 1:nt
        α = angles[t]
        Pmat = Matrix(generators[t])  # convert Pauli to matrix

        # Use in-place multiplication: W = U * P
        mul!(W, U, Pmat)

        # Update U in-place: U = cos(α/2)*U - i*sin(α/2)*W
        c, s = cos(α / 2), sin(α / 2)
        @inbounds @simd for i in eachindex(U)
            U[i] = c * U[i] - 1im * s * W[i]
        end
    end

    return U
end

function heisenberg(o::Pauli{N}; Jx, Jy, Jz, k) where N 
    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])
    # Loop over sites
    for ki in 1:k 
        for i in 1:N-1
            push!(generators, Pauli(N, X=[i, i + 1]))
            push!(parameters, Jx)
            push!(generators, Pauli(N, Y=[i, i + 1]))
            push!(parameters, Jy)
            push!(generators, Pauli(N, Z=[i, i + 1]))
            push!(parameters, Jz)
        end
    end

    return generators, parameters
end

function heisenberg_2D(o::Pauli{N}; Jx, Jy, Jz, k) where N 
    L = Int(sqrt(N))
    generators = Vector{Pauli{N}}()
    parameters = Vector{Float64}()

    index(i, j) = (mod1(i, L) - 1) * L + mod1(j, L)
    # this was for periodicboundary condition
    for ki in 1:k
        for i in 1:L, j in 1:L
            a = index(i, j)

            # Right neighbor (periodic in j)
            b = index(i, j + 1)
            push!(generators, Pauli(N, X=[a, b])); push!(parameters, Jx)
            push!(generators, Pauli(N, Y=[a, b])); push!(parameters, Jy)
            push!(generators, Pauli(N, Z=[a, b])); push!(parameters, Jz)

            # Down neighbor (periodic in i)
            b = index(i + 1, j)
            push!(generators, Pauli(N, X=[a, b])); push!(parameters, Jx)
            push!(generators, Pauli(N, Y=[a, b])); push!(parameters, Jy)
            push!(generators, Pauli(N, Z=[a, b])); push!(parameters, Jz)
        end
    end

    return generators, parameters
end

function heisenberg_2D_open(o::Pauli{N}; Jx, Jy, Jz, k) where N 
    L = Int(sqrt(N))
    generators = Vector{Pauli{N}}()
    parameters = Vector{Float64}()

    index(i, j) = (i - 1) * L + j  # no mod1, assume 1-based indexing

    for ki in 1:k
        for i in 1:L, j in 1:L
            a = index(i, j)

            # Right neighbor (open boundary)
            if j < L
                b = index(i, j + 1)
                push!(generators, Pauli(N, X=[a, b])); push!(parameters, Jx)
                push!(generators, Pauli(N, Y=[a, b])); push!(parameters, Jy)
                push!(generators, Pauli(N, Z=[a, b])); push!(parameters, Jz)
            end

            # Down neighbor (open boundary)
            if i < L
                b = index(i + 1, j)
                push!(generators, Pauli(N, X=[a, b])); push!(parameters, Jx)
                push!(generators, Pauli(N, Y=[a, b])); push!(parameters, Jy)
                push!(generators, Pauli(N, Z=[a, b])); push!(parameters, Jz)
            end
        end
    end

    return generators, parameters
end


function jw_transform(o::Pauli{N}, site) where N
    z_string = [i for i in 1:site-1]

    # p = PauliSum(N)
    p = Pauli(N, Z = z_string, X = [site]) + 1im * Pauli(N, Z=z_string, Y=[site])

    return 0.5*p
end

function fermi_hubbard_1D(o::Pauli{N}; t, U, k) where N
    Nsites = Int(N/2)
    generators = Vector{Pauli{N}}()
    parameters = Vector{Float64}()

    t_term = PauliSum(N)

    up(j) = 2*j - 1  
    dn(j) = 2*j      

    for ki in 1:k
        for j in 1:Nsites - 1
            # println("site: ", j)

            # α-spin c{i, α}†c{j, α} + h.c.
            i_a = jw_transform(o, up(j))
            j_a = jw_transform(o, up(j+1))
            
            t_term += i_a' * j_a
            t_term += j_a' * i_a

            # β-spin c{i, β}†c{j, β} + h.c. 
            i_b = jw_transform(o, dn(j))
            j_b = jw_transform(o, dn(j+1))
            t_term += i_b' * j_b
            t_term += j_b' * i_b 
            # println("Hopping")
            # display(temp)
        end

        for (pauli, coeff) in t_term.ops
            push!(generators, Pauli(pauli))
            push!(parameters, -t*coeff)
        end

        u_term = PauliSum(N)

        for j in 1:Nsites
            # interacting term
            i_a = jw_transform(o, up(j))
            i_b = jw_transform(o, dn(j))

            u_term += i_a'*i_a*i_b'*i_b
            # println("Interaction")
            # display(u_term)
        end

        for (pauli, coeff) in u_term.ops
            push!(generators, Pauli(pauli))
            push!(parameters, U*coeff)
        end
    end

    return generators, parameters
end 


function fermi_hubbard_2D(o::Pauli{N}; t, U, k) where N
    Nsites = Int(N/2)
    L = Int(sqrt(Nsites))
    generators = Vector{Pauli{N}}()
    parameters = Vector{Float64}()
    t_term = PauliSum(N)
    u_term = PauliSum(N)

    up(j) = 2*j - 1  
    dn(j) = 2*j
    linear_index(x, y) = (x-1)*L + y
    
    for ki in 1:k
        t_term = PauliSum(N)

        for x in 1:L
            for y in 1:L
                j = linear_index(x, y)
                if x < L
                    # down coupling
                    i = linear_index(x + 1, y)

                    # α-spin c{i, α}†c{j, α} + h.c.
                    i_a = jw_transform(o, up(j))
                    j_a = jw_transform(o, up(i))

                    t_term += i_a' * j_a
                    t_term += j_a' * i_a

                    # β-spin c{i, β}†c{j, β} + h.c. 
                    i_b = jw_transform(o, dn(j))
                    j_b = jw_transform(o, dn(i))
                    t_term += i_b' * j_b
                    t_term += j_b' * i_b  
                end  

                if y < L 
                    # α-spin c{i, α}†c{j, α} + h.c.
                    i = linear_index(x, y + 1)
                    # right side coupling
                    i_a = jw_transform(o, up(j))
                    j_a = jw_transform(o, up(i))

                    t_term += i_a' * j_a
                    t_term += j_a' * i_a

                    # β-spin c{i, β}†c{j, β} + h.c. 
                    i_b = jw_transform(o, dn(j))
                    j_b = jw_transform(o, dn(i))
                    t_term += i_b' * j_b
                    t_term += j_b' * i_b 
                end
            end
        end

        for (pauli, coeff) in t_term.ops
            push!(generators, Pauli(pauli))
            push!(parameters, -t*coeff)
        end 

        u_term = PauliSum(N)

        for j in 1:Nsites
            # interacting term
            i_a = jw_transform(o, up(j))
            i_b = jw_transform(o, dn(j))

            u_term += i_a'*i_a*i_b'*i_b
            # println("Interaction")
            # display(u_term)
        end

        for (pauli, coeff) in u_term.ops
            push!(generators, Pauli(pauli))
            push!(parameters, U*coeff)
        end
    end
    return generators, parameters, -t*t_term + U*u_term

end

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

# function central_spin(N, a, b, z, k)
#     generators = Vector{Pauli{N}}([])
#     parameters = Vector{Float64}([])

#     for ki in 1:k
#         for i in 2:N
#             push!(generators, Pauli(N, Z = [1],X = [i]))
#             push!(parameters, a)

#             push!(generators, Pauli(N, Z = [1], Z = [i]))
            
#             push!(generators, Pauli(N, Z = [i]))
#         end
#     end
# end

function central_spin(N, a, b, w)
    H = PauliSum(N)
    for i in 2:N
        H += a[i-1] * Pauli(N, Z = [1], X = [i])
        H += b[i-1] * Pauli(N, Z = [1, i])
        H += w * Pauli(N, Z = [i])
    end

    return H
end

function central_spin_gen(N, a, b, w, k)
    generators = Vector{Pauli{N}}([])
    parameters = Vector{Float64}([])

    for ki in 1:k
        for i in 2:N
            push!(generators, Pauli(N, Z = [1], X = [i]))
            push!(parameters,a[i-1])
            
            push!(generators, Pauli(N, Z = [1, i]))
            push!(parameters, b[i-1])

            push!(generators, Pauli(N, Z = [i]))
            push!(parameters, w)
        end
    end
    return generators, parameters
end
function local_folding(generators::Vector{Pauli{N}}, parameters, scaling_factor) where N
    # Adding U'U at rondom places in the generators

    places = rand(1:(length(generators)-1), scaling_factor)
    places = sort(places, rev = true)

    temp_gens = deepcopy(generators)
    temp_angles = deepcopy(parameters)

    for j in places
        add_pauli = temp_gens[j]
        add_angle = temp_angles[j]
        # temp_gens = vcat(temp_gens[1:j-1], add_pauli, add_pauli, temp_gens[j:end])
        # temp_angles = vcat(temp_angles[1:j-1], add_angle, -add_angle, temp_angles[j:end])

        insert!(temp_gens, j+1, -add_pauli)
        insert!(temp_gens, j+2, add_pauli)

        insert!(temp_angles, j+1, add_angle)
        insert!(temp_angles, j+2, add_angle)
      

    end
    return temp_gens, temp_angles
end
