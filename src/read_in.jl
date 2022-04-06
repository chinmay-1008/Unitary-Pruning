using UnitaryPruning, PyCall

np = pyimport("numpy")

ham_ops = Vector{PauliString{8}}()
ham_par = Vector{Float64}()
ansatz_ops = Vector{PauliString{8}}()
ansatz_par = Vector{Float64}()

for i in np.load("src/python/ham_ops.npy")
    push!(ham_ops, PauliString(i))
end
for i in np.load("src/python/ansatz_ops.npy")
    push!(ansatz_ops, PauliString(i))
end

for i in np.load("src/python/ansatz_par.npy")
    push!(ansatz_par, i)
end
for i in np.load("src/python/ham_par.npy")
    push!(ham_par, i)
end
