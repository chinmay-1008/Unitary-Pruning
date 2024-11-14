using UnitaryPruning
using Test
using BenchmarkTools 
using PauliOperators
using LinearAlgebra

# @testset "bfs_evolution" begin
function run_temp(N)
    ket = KetBitString(N, 0)

    op = [Pauli(N, X=[2,3], Y=[4], Z=[1,5]), Pauli(N, X=[3,4], Y=[1,2], Z=[5]), Pauli(N, X=[1], Y=[3,4], Z=[2,6]), Pauli(N, X=[1,2], Y=[6], Z=[3,4])]

    o = Pauli(N, X=[14,30,32], Y=[10,31], Z=[9,13,18,29,33])
    # i = 8
    # for o in op
    for i in 0:16
        α = i * π/32
        generators, parameters = UnitaryPruning.eagle_processor(o, α=α, k=8)
        # e = UnitaryPruning.bfs_evolution(generators, parameters, o, ket, thres=1e-4)
        e , nops, c_norm2 = UnitaryPruning.bfs_evolution(generators, parameters, PauliSum(o), ket, thresh=1e-3)

        println("Expectation value:             ", e)
        println("Norm 2 of the coefficients:    ", c_norm2)
        println()
    end
    return 0
end
energy = run_temp(127)
