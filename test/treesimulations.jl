rng = MersenneTwister(12)

@testset "tree branching" begin
    rng = MersenneTwister(12)
    @testset "initialisation" begin
        alivecells = SomaticEvolution.create_cells(TreeCell, WellMixed(), 10; rng)
        root = alivecells[1]
        @test root.data.mutations == 10
        alivecells = SomaticEvolution.create_cells(SimpleTreeCell, WellMixed(), 0; rng)
        root = alivecells[1]
        @test root.data.mutations == 0
        @test length(alivecells) == 1
        @test AbstractTrees.isroot(root)
        @test alivecells[1] == root
    end

    rng = MersenneTwister(12)
    @testset "run simulation" begin
        input = BranchingInput(
            Nmax=100, 
            birthrate=1, 
            deathrate=0.05,
            clonalmutations=0, 
            numclones=0,
            μ=10,
            mutationdist=:poisson
        )
        treemodule = runsimulation(TreeCell, input, rng).output
        treemodule_simple = runsimulation(SimpleTreeCell, input, rng).output
        root = getroot(SomaticEvolution.firstcellnode(treemodule))
        root_simple = getroot(SomaticEvolution.firstcellnode(treemodule_simple))
        @test length(SomaticEvolution.getalivecells(root)) == 100
        @test length(SomaticEvolution.getalivecells(root_simple)) == 100
        @test all(cellnode.data.alive for cellnode in treemodule.cells)
    end

    @testset "division" begin
        treemodule = SomaticEvolution.initialize(TreeCell, WellMixed, 0, 1)
        root = getsingleroot(treemodule.cells)
        _, nextID = SomaticEvolution.celldivision!(treemodule, 1, 1.1, 2, 10, :fixedtimedep, rng)
        @test nextID == 4
        @test length(treemodule.cells) == 2 #check number of alivecells equals pop size
        @test !(root in treemodule.cells) #check divided cell has been removed from alivecells
        @test root.data.mutations == 11 #check divided cell has correct number of mutations
        #check new cells have correct properties and are children of root
        for (i, cellnode) in enumerate(treemodule.cells)
            @test ischild(cellnode, root)
            @test cellnode.data.alive
            @test cellnode.data.birthtime ≈ 1.1 atol=0.001
            @test cellnode.data.mutations == 0
            @test cellnode.data.clonetype == root.data.clonetype
            @test cellnode.data.id == i + 1
        end

        #asymmetric division
        _, nextID = SomaticEvolution.celldivision!(treemodule, 2, 1.1, 3, 10, :fixedtimedep, rng; nchildcells=1)
        @test length(treemodule.cells) == 2 #check number of alivecells equals pop size
        @test !(root in treemodule.cells) #check divided cell has been removed from alivecells
        @test treemodule.cells[end] == root.right.left
        @test isnothing(root.right.right) 
        @test treemodule.cells[end].data.id == 3


    end


    @testset "death" begin
        treemodule = SomaticEvolution.initialize(TreeCell, WellMixed, 0, 1)
        alivecells = treemodule.cells
        root = getsingleroot(alivecells)
        _, nextID = SomaticEvolution.celldivision!(treemodule, 1, 1.1, 2, 10, :fixedtimedep, rng)
        SomaticEvolution.celldeath!(treemodule, 1, 1.5, 10, :fixedtimedep, rng)
        @test length(alivecells) == 1
        @test !(root.left.left.data.alive) #check dead cell has correct properties
        @test root.left.left.data.birthtime ≈ 1.5 atol=0.0001
        @test root.left.left.data.mutations == 0
        @test root.left.left.data.clonetype == root.left.data.clonetype
        @test root.left.data.mutations == 4 #check dead cell has correct number mutations
        @test isnothing(root.left.right) #check there is no right child of dead cell

    end

    @testset "death simple" begin
        treemodule = SomaticEvolution.initialize(SimpleTreeCell, WellMixed, 0, 1)
        alivecells = treemodule.cells
        root = getsingleroot(alivecells)
        _, nextID = SomaticEvolution.celldivision!(treemodule, 1, 1.1, 2, 10, :fixedtimedep, rng)
        SomaticEvolution.celldeath!(treemodule, 1, 1.5, 10, :fixedtimedep, rng)
        @test length(alivecells) == 1
        _, nextID = SomaticEvolution.celldivision!(treemodule, 1, 1.1, 2, 10, :fixedtimedep, rng)
        idx = 2
        deadcellnode = alivecells[idx]
        SomaticEvolution.celldeath!(treemodule, idx)
        @test !AbstractTrees.intree(deadcellnode, root) #check dead cell is not in tree
        @test !(deadcellnode in alivecells) #check dead cell is not in alivecells vector

    end

end

@testset "tree moran" begin
    input = MoranInput(N=4, tmax=5, mutationdist=:fixed)
    treemodule = runsimulation(SimpleTreeCell, input, rng).output
    @test length(treemodule.cells) == length(treemodule) == 4
    @test popsize(getroot(treemodule.cells)) == 4
    treemodule = runsimulation(TreeCell, input, rng).output
    @test length(treemodule.cells) == length(treemodule) == 4
    @test popsize(getroot(treemodule.cells)) == 4
end


#create a small tree structure with 5 nodes and 3 alive cells
# -root -|- leftcell -|- leftleftcell
#        |            |- leftrightcell
#        |- rightcell

function make_tree(::Type{T}) where T <: SomaticEvolution.AbstractTreeCell
    rootcell = T(id=1, birthtime=0.0, mutations=0)
    leftcell = T(id=2, birthtime=1.5355542835848743, mutations=5)
    rightcell = T(id=3, birthtime=1.5355542835848743, mutations=10)
    leftleftcell = T(id=4, birthtime=1.6919799338516708, mutations=13)
    leftrightcell = T(id=5, birthtime=1.6919799338516708, mutations=17)
    root = BinaryNode(rootcell)
    leftchild!(root, leftcell)
    rightchild!(root, rightcell)
    leftchild!(root.left, leftleftcell)
    rightchild!(root.left, leftrightcell)
    return root
end
tree_root = make_tree(TreeCell)
tree_root_simple = make_tree(SimpleTreeCell)

@testset "basic tree" begin
    for root in [tree_root, tree_root_simple]
        @test AbstractTrees.nextsibling(root.left.left) == root.left.right
        @test isnothing(AbstractTrees.nextsibling(root.left.right))    
        @test AbstractTrees.prevsibling(root.left.right) == root.left.left  
        @test isnothing(AbstractTrees.prevsibling(root.left.left))    
        @test isnothing(AbstractTrees.nextsibling(root))
    end
end

@testset "prune tree" begin
    #check that prune_tree! works correctly, i.e. if it removes a node, leaving a parent
    #node with no children, it removes that parent node etc.

    root = make_tree(TreeCell)
    alivecells = [root.left.left, root.left.right, root.right]
    SomaticEvolution.prune_tree!(root.left.left)
    @test !intree(root.left.left, root)
    SomaticEvolution.prune_tree!(root.left.right)
    @test !intree(root.left, root)

end

@testset "tree statistics" begin
    for root in [tree_root, tree_root_simple]
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
end

#create a second tree with 3 nodes and 1 alive cell
function make_tree2(::Type{T}) where T <: SomaticEvolution.AbstractTreeCell
    root = BinaryNode(T(id=1, birthtime=2.0, mutations=0))
    leftchild!(root, T(id=1, birthtime=3.0, mutations=5))
    rightchild!(root, T(id=1, birthtime=3.0, mutations=11))
    return root
end

tree_root2 = make_tree2(TreeCell)
tree_root_simple2 = make_tree2(SimpleTreeCell)
tree_roots = vcat(tree_root, tree_root2)
tree_roots_simple = vcat(tree_root_simple, tree_root_simple2)
alivecells_2roots = [tree_root.right, tree_root.left.left, tree_root.left.right, tree_root2.left, tree_root2.right]
alivecells_2roots_simple = [tree_root_simple.right, tree_root_simple.left.left, tree_root_simple.left.right, tree_root_simple2.left, tree_root_simple2.right]


@testset "multiple roots" begin
    for (root, root2, roots, alive2roots) in [
        (tree_root, tree_root2, tree_roots, alivecells_2roots), 
        (tree_root_simple, tree_root_simple2, tree_roots_simple, alivecells_2roots_simple)
    ]
        @test getroot([root2.left, root2.right]) == [root2]
        @test Set(getroot(alive2roots)) == Set(roots)
        @test Set(getalivecells(roots)) == Set(alive2roots)
        @test isnothing(getsingleroot(alive2roots)) 
        @test getsingleroot([root2.left, root2.right]) == root2
    end
end

@testset "MRCA" begin
    root = make_tree(SimpleTreeCell)
    alivecells = collect(Leaves(root))
    @test findMRCA(alivecells) == root
    @test findMRCA(root.left.left, root.left.right) == root.left
    @test findMRCA(root.left.left, root.right) == root
end