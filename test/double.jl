using PauliOperators
using LinearAlgebra
using Random
using Plots
using Printf

Random.seed!(1234)


function newclip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(p->(abs(p.second) ≥ thresh) , ps)
end


function isantihermitian(p::PauliSum{N,T}) where {N,T}
    isanti = true
    for coeff in values(p)
        # Check if coefficient is (approximately) purely imaginary
        isanti = isanti && isapprox(real(coeff), 0, atol=1e-16)
    end
    return isanti
end


function heisenberg_1d(N, jx, jy, jz)
    hammy = PauliSum(N)
    for i in 1:N-1  
        hammy += jx * Pauli(N, X=[i, i+1])
        hammy += jy * Pauli(N, Y=[i, i+1])
        hammy += jz * Pauli(N, Z=[i, i+1])
    end
    return hammy
end


function ising_1D(N, jx, jz)
    hammy = PauliSum(N)

    for i in 1:N-1
        hammy += -jz * Pauli(N, Z = [i, i+1])
    end

    for i in 1:N
        hammy += -jx * Pauli(N, X=[i])  
    end
    return hammy
end


function evolve(P::PauliSum{N,T}, G::Pauli{N}, dt) where {N,T}
    _cos = cos(dt*coeff(G))
    _sin = 1im*sin(dt*coeff(G))
    out = deepcopy(P) 
    sin_branch = PauliSum(N)
    for (p,c) in P
        if PauliOperators.commute(p,PauliBasis(G)) == false
            out[p] *= _cos
            sum!(sin_branch, c*_sin*p*PauliBasis(G))
        end
    end
    sum!(out, sin_branch)
    return out
end


function evolve_full(P::Union{PauliSum{N, T}, Pauli{N}}, G::PauliSum{N, T}, dt) where {N, T}
# We want to do  ∑_j [∏e^(iθG)] c_jP_j [∏e^(-iθG)]
    out =  PauliSum(N)
    generators = []
    params = []

    for (g, c) in G
        # put a sanity check
        if isapprox(real(c), 0; atol=1e-16)
            # display("real coefficient found: $c")
            push!(generators, g)
            push!(params, dt*imag(c))
        end

        # push!(generators, g)
        # push!(params, c)
    end
    nt = length(generators)

    for (p, coeff) in P
        o_transformed = coeff*PauliSum(p) 
        # sin_branch = PauliSum(N)
        for t in 1:nt
            sin_branch = PauliSum(N)
            for (o, o_coeff) in o_transformed
                g = generators[t]
                if PauliOperators.commute(o, g) == false
                    o_transformed[o] = o_coeff * cos(2 * params[t])
                    sum!(sin_branch, 1im * o_coeff * g * o * sin(2 * params[t]) )
                    # o_transformed[o] = o_coeff * cos(params[t])
                    # sum!(sin_branch, 1im * o_coeff * g * o * sin(params[t]) )
                end
            end
            sum!(o_transformed, sin_branch)
            # newclip!(o_transformed, thresh=1e-11)
        end
        sum!(out, o_transformed)
    end
    return out
end 


function evolve_full_new(P::PauliSum{N, T}, G::PauliSum{N, T}, dt) where {N, T}
    # We want to do  ∑_j [∏e^(iθG)] c_jP_j [∏e^(-iθG)]
        out =  PauliSum(N)
        generators = []
        params = []
    
        for (g, c) in G
            # put a sanity check
            if isapprox(real(c), 0; atol=1e-13)
                # display("real coefficient found: $c")
                push!(generators, g)
                push!(params, dt*imag(c))
            else
                display(c)
                throw(ErrorException("Not antihermitian!"))
            end
        end
        # display(params)
        nt = length(generators)
        
        o_transformed = deepcopy(P) 
    
        for t in 1:nt
            g = generators[t]
            sin_branch = PauliSum(N)
            for (o, o_coeff) in o_transformed
                if PauliOperators.commute(o, g) == false
                    o_transformed[o] = o_coeff * cos(2 * params[t])
                    sum!(sin_branch, 1im * o_coeff * g * o * sin(2 * params[t]) )
                end
            end
            sum!(o_transformed, sin_branch)
        end
        return o_transformed
        # sum!(out, o_transformed)
        # return out
    end 
    

function largest_term(ps::PauliSum)
    best = nothing
    maxval = 0.0
    coeff = 0.0
    for (p, c) in ps
        v = abs(c)
        if v > maxval
            maxval = v
            best = p
            coeff = c
        end
    end
    return coeff*PauliSum(best)
end

function diagonal_paulisum(N::Int)
    op = PauliSum(N) 

    for mask in 1:(2^N - 1)
        zsites = [i for i in 1:N if (mask >> (i-1)) & 1 == 1]
        op += randn()*Pauli(N, Z=zsites)
    end
    return op
end


function unitary_matrix(G::PauliSum{N}, dt) where {N}
    # Build Hermitian generator: angles are real
    θs = [dt * imag(c) for (_, c) in G]
    gs = [g for (g, _) in G]

    U = I # start with identity
    for (θ, g) in zip(θs, gs)
        # exp(i θ g) = cos(θ)I + i sin(θ) g
        Gmat = Matrix(g)        # convert Pauli operator to dense matrix
        U = (cos(θ) * I + 1im * sin(θ) * Gmat) * U
    end
    return U
end

# ```
# We want to do the double bracket evolution
# dH/dt = -[H, [D, H]] where D is the diagonal terms of the Hamiltonian.
# THis gives rise to H(t + dt) = U H(t) U', where U = e^(dt*[D, H]), and [D, H] is anti-hermitian
# In small limit of dt we can write this equation as a similarity transformation or heisenberg evolution
# ```

function test_evolve()
    N = 2
    H = heisenberg_1d(N, 0.8, 0.9, 0.9)
    o = PauliSum(Pauli(N, Z = [1]))
    dt = 1

    k = 10
    for i in 1:k
        o = evolve_full(o, H, dt)
    end
    display(o)
    println("NEW")
    display(expectation_value(o, Ket(N, 0)))

          
    ket = Ket(N, 0) 
    o = Pauli(N, Z=[1])
    println("BFS")
    generators, parameters = UnitaryPruning.heisenberg(o, Jx = 0.8, Jy = 0.9,Jz = 0.9, k=k)
    ei , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=1e-10)       
    display(ei)

    println("EXACT")
    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    o_mat = Matrix(o)
    m = diag(U'*o_mat*U)
    display(m[1])
end





function run()
    N = 2
    H = rand(PauliSum{N}, n_paulis = 50)
    # H = heisenberg_1d(N, 1, 1, 5)
    # H = ising_1D(N, 1, 2)
    H += H'

    eigval, _ = eigen(Matrix(H))
    t_steps = 10000
    dt = 1e-4

    # D = diagonal_paulisum(N)
    # display(D)
    # return D
    errs = []
    times = Float64[]

    for i in 1:t_steps

        D = PauliSum(N)
        for (p, c) in H
            if p.x == 0
                sum!(D, c*p)
            end
        end
        com = D * H - H * D

        # newclip!(com, thresh = 1e-6)

        # u = unitary_matrix(com, dt)
        # display(norm(u' * u - I))

        # com = largest_term(com)

        H = evolve_full(H, com, dt)

        # newclip!(com, thresh = 1e-10)
        if i % 100 == 0
            # display(H)
            mat = Matrix(H)
            err = eigval - sort(diag(mat), by= real)
            # err = norm(mast - Diagonal(diag(mat)))

            # display(mat)
            ne = norm(err)
            println("Time step:", i,"  ", ne)

            push!(errs, ne)
            push!(times, i * dt)   
        end
    end

    println("Eigenvalues")
    display(eigval)
    println("Evolved Hamiltoian")
    display(sort(real(diag(Matrix(H)))))
    # display(Matrix(H))

    eigval, _ = eigen(Matrix(H))
    println("Evolved Hamiltoian Diagonalization")
    display(eigval)

    plot(times, errs, xlabel="time", ylabel="‖error‖", lw=2, legend=false, markers=false, dpi=300,
    title = "Norm of error of diagonal terms; dt = $dt")
    savefig("test/error_single.png")

    return
end

function run_matrix()

    N = 3
    H = rand(PauliSum{N}, n_paulis = 50)
    # H = heisenberg_1d(N, 1, 1, 5)
    # H = ising_1D(N, 1, 2)
    H += H'

    H = Matrix(H)

    eigval, _ = eigen(H)
    t_steps = 100000
    dt = 1e-5
    # display(eigval)
    # return
    errs = []
    times = Float64[]
    # D = diagonal_paulisum(N)
    # D = Matrix(D)

    for i in 1:t_steps
        # Step 1: Extract diagonal part of H
        D = Diagonal(diag(H))
    
        # Step 2: Build generator G = [D, H]
        G = D*H - H*D   # commutator
    
        # Step 3: Exponentiate generator
        U = exp(dt * G)
        # display(U)
        
    
        # Step 4: Update H by similarity transform
        H = U * H * U'
        # display(U')
        if i % 1000 == 0
            # ne = norm(sort(abs.(eigval)) - sort(abs.(diag(H))))
            ne = norm(H - Diagonal(diag(H)))

            println(ne)
            push!(errs, ne)
            push!(times, i * dt)     
        end

    end
    println("Eigenvalues")
    display(eigval)
    println("Evolved Hamiltoian")
    display(sort(real(diag(H))))
    eigval, _ = eigen(H)
    println("Evolved Hamiltoian Diagonalization")
    display(eigval)

    plot(times, errs, xlabel="time", ylabel="‖error‖", lw=2, legend=false, markers=false, dpi=300,
    title = "Norm of error of diagonal terms; dt = $dt")
    savefig("test/error_full_matrix.png")
end

function run_discrete()
    N = 3
    H = rand(PauliSum{N}, n_paulis = 50)
    # H = heisenberg_1d(N, 1, 1, 5)
    # H = ising_1D(N, 1, 2)
    H += H'

    # H = Matrix(H)

    eigval, _ = eigen(Matrix(H))
    t_steps = 100000
    dt = 1e-6
    # display(eigval)
    # return
    errs = []
    times = Float64[]
    # D = diagonal_paulisum(N)
    # D = Matrix(D)

    for i in 1:t_steps
        # Step 1: Extract diagonal part of H
        D = PauliSum(N)
        for (p, c) in H
            if p.x == 0
                sum!(D, c*p)
            end
        end
        com = D * H - H * D
    
        # Step 3: create the double commutator
        dcom = H * com - com * H
        # display(U)
        
    
        # Step 4: Update H by similarity transform
        H = H - dt * dcom
        # display(U')
        if i % 1000 == 0
            # ne = norm(sort(abs.(eigval)) - sort(abs.(diag(H))))
            mat = Matrix(H)
            ne = norm(mat - Diagonal(diag(mat)))

            println("Time step:", i,"  ", ne)
            push!(errs, ne)
            push!(times, i * dt)     
        end

    end
    println("Eigenvalues")
    display(eigval)
    println("Evolved Hamiltoian")
    display(sort(real(diag(Matrix(H)))))
    eigval, _ = eigen(Matrix(H))
    println("Evolved Hamiltoian Diagonalization")
    display(eigval)

    plot(times, errs, xlabel="time", ylabel="‖error‖", lw=2, legend=false, markers=false, dpi=300,
    title = "Norm of error of diagonal terms; dt = $dt")
    savefig("test/error_full_discrete.png")
    
end

function compare_flows()

    N = 3
    t_steps = 10000
    dt = 1e-3

    # Initial Hermitian Hamiltonian
    H0 = rand(PauliSum{N}, n_paulis = 50)
    H0 += H0'
    eigval, _ = eigen(Matrix(H0))

    errs_pauli = Float64[]
    errs_matrix = Float64[]
    errs_discrete = Float64[]
    times = Float64[]

    # --- Method 1: Pauli propagation ---
    H_pauli = copy(H0)

    # --- Method 2: Dense matrix similarity ---
    H_matrix = Matrix(H0)
    # D_fixed = Matrix(diagonal_paulisum(N))  # fixed diagonal operator

    # --- Method 3: Discrete double-commutator ---
    H_disc = copy(H0)

    for i in 1:t_steps
        # ----- Method 1 -----
        Dp = PauliSum(N)
        for (p, c) in H_pauli
            if p.x == 0
                sum!(Dp, c*p)
            end
        end
        com = Dp * H_pauli - H_pauli * Dp
        # com = largest_term(com)        # pruning to single Pauli generator
        H_pauli = evolve_full_new(H_pauli, com, dt)

        # ----- Method 2 -----
        Dm = Diagonal(diag(H_matrix))
        G = Dm * H_matrix - H_matrix * Dm
        U = exp(dt * G)
        H_matrix = U * H_matrix * U'

        # ----- Method 3 -----
        com_d = Dp * H_disc - H_disc * Dp
        dcom = H_disc * com_d - com_d * H_disc
        H_disc = H_disc - dt * dcom

        if i % 100 == 0
            push!(times, i * dt)

            mat1 = Matrix(H_pauli)
            mat2 = H_matrix
            mat3 = Matrix(H_disc)

            # push!(errs_pauli,  norm(mat1 - Diagonal(diag(mat1))))
            # push!(errs_matrix, norm(mat2 - Diagonal(diag(mat2))))
            # push!(errs_discrete, norm(mat3 - Diagonal(diag(mat3))))

            push!(errs_pauli,  norm(eigval - sort(diag(mat1), by = real)))
            push!(errs_matrix, norm(eigval - sort(diag(mat2), by = real)))
            push!(errs_discrete, norm(eigval - sort(diag(mat3), by = real)))

            println("Time step:", i)

        end
    end

    # Plot all on one figure
    plot(times, errs_pauli, lw=2, label="Pauli propagation", dpi = 300, yscale=:log10,
    )
    plot!(times, errs_matrix, lw=2, label="Matrix formalism")
    plot!(times, errs_discrete, lw=2, label="Discrete double commutator")
    xlabel!("time")
    ylabel!("‖error‖")
    title!("Double-bracket flow comparison (N=$N, dt=$dt)")
    savefig("test/compare_flows_eigs.png")

    # println("\nInitial spectrum:")
    # init_spec = sort(eigval)
    # for i in eachindex(init_spec)
    #     @printf("%20.12f\n", init_spec[i])
    # end    
    println("Initial spectrum vs Final spectrum vs Final diagonals")

    methods = [
        ("Pauli Propagation", Matrix(H_pauli)),
        ("Matrix Formalism", H_matrix),
        ("Discretised ODE", Matrix(H_disc))
    ]
    
    for (name, Hf) in methods
        println("\n$name")
        println(@sprintf("%20s | %20s | %20s", "eigvals(H(0))", "eigvals(H(T))", "diag(H(T))"))
        println("-"^70)
    
        init_spec = sort(eigval)
        final_spec = sort(real(eigen(Hf).values))
        diagvals  = sort(real(diag(Hf)))
    
        for i in eachindex(init_spec)
            @printf("%20.12f | %20.12f | %20.12f\n", init_spec[i], final_spec[i], diagvals[i])
        end
    end
    
    return
end


compare_flows()
# @time run()
# @time run_matrix()
# @time run_discrete()
# test_evolve()