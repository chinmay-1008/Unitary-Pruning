using PauliOperators
using LinearAlgebra
using Random
using Plots

Random.seed!(1234)

# clipping small coefficients
function newclip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(pair -> abs(pair[2]) ≥ thresh, ps)
end

# evolve with PauliSum generator
function evolve_full(P::Union{PauliSum{N, T}, Pauli{N}}, G::PauliSum{N, T}, dt) where {N, T}
    Psum = isa(P, Pauli) ? PauliSum(P) : deepcopy(P)
    out = PauliSum(N)

    for (p, coeff) in Psum
        o_transformed = PauliSum(p)
        for (d, d_coeff) in G
            next_op = PauliSum(N)
            for (o, o_coeff) in o_transformed
                if PauliOperators.commute(o, d)
                    sum!(next_op, o_coeff * o)
                else
                    cos_part = o_coeff * cos(dt * d_coeff)
                    sin_part_coeff = o_coeff * (1im * sin(dt * d_coeff))
                    sum!(next_op, cos_part * o)
                    sum!(next_op, sin_part_coeff * (d * o))
                end
            end
            o_transformed = next_op
            newclip!(o_transformed; thresh=1e-12)
        end
        newclip!(o_transformed; thresh=1e-12)
        sum!(out, o_transformed * coeff)
    end
    newclip!(out; thresh=1e-12)
    return out
end

# extract diagonal-only Pauli terms
function diagonal_part(H::PauliSum{N}) where {N}
    D = PauliSum(N)
    for (p, c) in H
        if p.x == 0
            sum!(D, c * p)
        end
    end
    return D
end

function run_animation(; N=3, steps=200, dt=1e-3)
    # start Hamiltonian
    H = rand(PauliSum{N}, n_paulis=20)
    H += H'
    evals = eigvals(Matrix(H))

    println("Initial eigenvalues: ", sort(real(evals)))

    # collect frames
    frames = @gif for i in 1:steps "test/pauli_propagation.gif"
        # diagonal projector
        D = diagonal_part(H)
        com = H * D - D * H
        H = evolve_full(H, com, dt)

        matH = real(Matrix(H))

        heatmap(matH,
            c=:RdGy, #clim=(-maximum(abs.(matH)), maximum(abs.(matH))),
            title="Step $i", xlabel="", ylabel="", colorbar=false, aspect_ratio=1)
    end  # record every 2nd step to keep animation small

    return frames
end

anim = run_animation(N=4, steps=150, dt=1e-3)
