alivecells, root, N, nextID = initialize_branching_tree(0, 1)
@testset "tree initialisation" begin
    @test length(alivecells) == 1 == N
    @test AbstractTrees.isroot(root)
    @test nextID == root.data.id + 1
    @test alivecells[1] == root
end

input = BranchingInput(
    Nmax=100, 
    b=1, 
    d=0.0,
    clonalmutations=0, 
    numclones=0,
    μ=1,
    mutationdist=:fixed
)
rng = MersenneTwister(12)
alivecells, root = run1simulation_tree(input, rng)
@testset "tree basic" begin
    
end

tree = BinaryNode(SimpleCell(1, true, 0.0, 0, 1))
leftchild(SimpleCell(2, true, 1.5355542835848743, 5, 1), tree)
rightchild(SimpleCell(3, true, 1.5355542835848743, 10, 1), tree)
leftchild(SimpleCell(4, true, 1.6919799338516708, 13, 1), tree.left)
rightchild(SimpleCell(5, true, 1.6919799338516708, 17, 1), tree.left)

@testset "tree statistics" begin
    @test mutations_per_cell(tree, includeclonal=true) == [18, 22, 10] 
    @test celllifetimes(tree, excludeliving=true) ≈ [1.53555, 0.15643] atol=0.01
    @test celllifetimes(tree, excludeliving=false) ≈ [1.53555, 0.15643, 0.0, 0.0, 0.15643] atol=0.01
    @test time_to_MRCA(tree.left.left, tree.right, 2.0) ≈ 2.0 - 1.5355542835848743
    @test time_to_MRCA(tree.left.left, tree.left.right, 2.0) ≈ 2.0 - 1.6919799338516708
    @test coalescence_times(tree) ≈ [0.0, 1.6919799338516708 - 1.5355542835848743, 1.6919799338516708 - 1.5355542835848743]
    @test pairwisedistance(tree.left.left, tree.right) == 28
    @test pairwisedistance(tree.left.left, tree.left.right) == 30
    @test pairwise_fixed_differences(tree) == Dict(28=>1, 30=>1, 32=>1)

end