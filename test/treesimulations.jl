rng = MersenneTwister(12)

@testset "tree branching" begin
    @testset "initialisation" begin
        alivecells = initialize_tree(10)
        root = getroot(alivecells)
        @test root.data.mutations == 10
        alivecells = initialize_tree(0)
        root = getroot(alivecells)
        @test root.data.mutations == 0
        @test length(alivecells) == 1
        @test AbstractTrees.isroot(root)
        @test alivecells[1] == root
    end

    @testset "run simulation" begin
        alivecells = initialize_tree(0)
        root = getroot(alivecells)
        input = BranchingInput(
            Nmax=100, 
            b=1, 
            d=0.05,
            clonalmutations=0, 
            numclones=0,
            μ=10,
            mutationdist=:poisson
        )
        alivecells = run1simulation_tree(input, rng)
        root = getroot(alivecells)
        @test length(SomaticEvolution.getalivecells(root)) == 100
        @test all(cellnode.data.alive for cellnode in alivecells)
    end

    @testset "division" begin
        alivecells, root = initialize_tree(0)
        N, nextID = SomaticEvolution.celldivision!(alivecells, 1, 1, 1.1, 2, 10, :fixedtimedep, rng)
        @test N == 2
        @test nextID == 4
        @test length(alivecells) == 2 #check number of alivecells equals pop size
        @test !(root in alivecells) #check divided cell has been removed from alivecells
        @test root.data.mutations == 11 #check divided cell has correct number of mutations
        #check new cells have correct properties and are children of root
        for (i, cellnode) in enumerate(alivecells)
            @test ischild(cellnode, root)
            @test cellnode.data.alive
            @test cellnode.data.birthtime ≈ 1.1 atol=0.001
            @test cellnode.data.mutations == 0
            @test cellnode.data.clonetype == root.data.clonetype
            @test cellnode.data.id == i + 1
        end
    end

    @testset "death" begin
        alivecells, root = initialize_tree(0)
        N, nextID = SomaticEvolution.celldivision!(alivecells, 1, 1, 1.1, 2, 10, :fixedtimedep, rng)
        N = SomaticEvolution.celldeath!(alivecells, 1, 2, 1.5, 10, :fixedtimedep, rng)
        @test N == 1
        @test !(root.left.left.data.alive) #check dead cell has correct properties
        @test root.left.left.data.birthtime ≈ 1.5 atol=0.0001
        @test root.left.left.data.mutations == 0
        @test root.left.left.data.clonetype == root.left.data.clonetype
        @test root.left.data.mutations == 4 #check dead cell has correct number mutations
        @test !(isdefined(root.left, :right)) #check there is no right child of dead cell

    end

end


#create a small tree structure with 5 nodes and 3 alive cells
rootcell = SimpleCell(1, true, 0.0, 0, 1)
leftcell = SimpleCell(2, true, 1.5355542835848743, 5, 1)
rightcell = SimpleCell(3, true, 1.5355542835848743, 10, 1)
leftleftcell = SimpleCell(4, true, 1.6919799338516708, 13, 1)
leftrightcell = SimpleCell(5, true, 1.6919799338516708, 17, 1)
root = BinaryNode(rootcell)
leftchild(leftcell, root)
rightchild(rightcell, root)
leftchild(leftleftcell, root.left)
rightchild(leftrightcell, root.left)

@testset "tree statistics" begin
    @test mutations_per_cell(root, includeclonal=true) == [18, 22, 10] 
    @test celllifetimes(root, excludeliving=true) ≈ [1.53555, 0.15643] atol=0.01
    @test celllifetimes(root, excludeliving=false) ≈ [1.53555, 0.15643, 0.0, 0.0, 0.15643] atol=0.01
    @test time_to_MRCA(root.left.left, root.right, 2.0) ≈ 2.0 - 1.5355542835848743
    @test time_to_MRCA(root.left.left, root.left.right, 2.0) ≈ 2.0 - 1.6919799338516708
    @test coalescence_times(root) ≈ [0.0, 1.6919799338516708 - 1.5355542835848743, 1.6919799338516708 - 1.5355542835848743]
    @test pairwisedistance(root.left.left, root.right) == 28
    @test pairwisedistance(root.left.left, root.left.right) == 30
    @test pairwise_fixed_differences(root) == Dict(28=>1, 30=>1, 32=>1)
    @test Set(getalivecells(root)) == Set([root.right, root.left.left, root.left.right])
    @test endtime(root.left.left) === nothing
    @test endtime(root) ≈ 1.5355542835848743 atol=1e-6
end