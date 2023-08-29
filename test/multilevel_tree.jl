@testset "simulate to fixed size" begin
    birthrate = 0.1
    deathrate = 0.0
    moranrate = 0.01
    asymmetricrate = 0.01
    modulesize = 4
    branchinitsize = 1
    branchrate = 0.01
    tmax = 1000
    maxmodules = 5
    mutationdist = :fixed
    μ = 1
    modulesplitting_replacement=false

    input = MultilevelBranchingInput(;
        birthrate, 
        deathrate, 
        moranrate, 
        asymmetricrate,
        modulesize, 
        branchinitsize, 
        branchrate, 
        tmax,
        maxmodules, 
        μ, 
        mutationdist
    )

    rng = MersenneTwister(1)
    moranincludeself=true
    population = SomaticEvolution.initialize_population(TreeCell, WellMixed, input.clonalmutations, SomaticEvolution.getNinit(input); rng)
    nextID = SomaticEvolution.getnextID(population) #get id of next cell
    nextmoduleID = 2
    @test nextID == 2
    t = 0
    rng = MersenneTwister(1)

    #check first update is correct
    transitionrates = SomaticEvolution.get_transitionrates(
        population, 
        birthrate, 
        deathrate, 
        moranrate, 
        asymmetricrate,
        branchrate, 
        modulesize
    ) 
    #only tranisition with non-zero rate is birth which has rate birthrate=0.1
    @test transitionrates == [0.0, 0.0, 0.1, 0.0, 0.0] 
    population, t, nextID = 
        SomaticEvolution.update_population!(
            population, 
            birthrate, 
            deathrate, 
            moranrate, 
            asymmetricrate,
            branchrate, 
            modulesize, 
            branchinitsize,
            modulesplitting_replacement,
            t,
            nextID, 
            nextmoduleID,
            μ,
            mutationdist, 
            tmax,
            maxmodules,
            moranincludeself,
            rng
    )
    #popsize is now 2
    @test length(population[1]) == length(population[1].cells) == 2
    #each cell has one mutation
    @test population[1].cells[1].data.mutations == population[1].cells[2].data.mutations == 1

    #check runsimulation function
    rng = MersenneTwister(12)
    input = MultilevelBranchingInput(
            modulesize=4, 
            fixedmu=true, 
            birthrate=0.1, 
            deathrate=0.01,
            moranrate=0.01, 
            asymmetricrate=0.01, 
            clonalmutations=0, 
            tmax=20*365, 
            maxmodules=10,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = runsimulation(SimpleTreeCell, input, rng)
    @test length(population) == 10
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    population = runsimulation(TreeCell, input, rng)
    @test length(population) == 10
    input = MultilevelBranchingInput(
            modulesize=4, 
            fixedmu=true, 
            birthrate=1, 
            deathrate=0,
            moranrate=0.1, 
            clonalmutations=0, 
            tmax=365*100, 
            maxmodules=10,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = runsimulation(SimpleTreeCell, input, rng)
    @test age(population) <= 365*50
    population = runsimulation(TreeCell, input, rng)
end

# 1 -- 2 -- 3 -- 5 -- 7 -- 8 -- 16 -- 19
#      |                         | -- 20 -- 24
#      |                               | -- 25
#      |
#      | -- 4 -- 6 -- 9 -- 10 -- 12 -- 14 -- 17 -- 21 -- 26
#                     |     |                       | -- 27
#                     |     |
#                     |     | -- 13 -- 28 
#                     |
#                     | -- 11 -- 15 -- 18 -- 22 -- 29
#                                       |     | -- 30
#                                       |         
#                                       | -- 23

#region Define tree structured population
cells = [SimpleTreeCell(i, i*0.1, i, 1) for i in 1:30]
root = BinaryNode(cells[1])
leftchild!(root, cells[2])
leftchild!(root.left, cells[3])
leftchild!(root.left.left, cells[5])
leftchild!(root.left.left.left, cells[7])
leftchild!(root.left.left.left.left, cells[8])
mrca1 = leftchild!(root.left.left.left.left.left, cells[16])
node19 = leftchild!(root.left.left.left.left.left.left, cells[19])
rightchild!(root.left.left.left.left.left.left, cells[20])
node24 = leftchild!(root.left.left.left.left.left.left.right, cells[24])
node25 = rightchild!(root.left.left.left.left.left.left.right, cells[25])

rightchild!(root.left, cells[4])
leftchild!(root.left.right, cells[6])
leftchild!(root.left.right.left, cells[9])
mrca2 = leftchild!(root.left.right.left.left, cells[10])
leftchild!(root.left.right.left.left.left, cells[12])
leftchild!(root.left.right.left.left.left.left, cells[14])
leftchild!(root.left.right.left.left.left.left.left, cells[17])
leftchild!(root.left.right.left.left.left.left.left.left, cells[21])
node26 = leftchild!(root.left.right.left.left.left.left.left.left.left, cells[26])
node27 = rightchild!(root.left.right.left.left.left.left.left.left.left, cells[27])
rightchild!(root.left.right.left.left.left, cells[13])
node28 = leftchild!(root.left.right.left.left.left.right, cells[28])
rightchild!(root.left.right.left.left, cells[11])
leftchild!(root.left.right.left.left.right, cells[15])
mrca3 = leftchild!(root.left.right.left.left.right.left, cells[18])
leftchild!(root.left.right.left.left.right.left.left, cells[22])
node29 = leftchild!(root.left.right.left.left.right.left.left.left, cells[29])
node30 = rightchild!(root.left.right.left.left.right.left.left.left, cells[30])
node23 = rightchild!(root.left.right.left.left.right.left.left, cells[23])

structure = WellMixed()
module1 = SomaticEvolution.TreeModule(Union{BinaryNode{SimpleTreeCell}, Nothing}[node19, node24, node25], 5.0, [0.0, 2.0], SomaticEvolution.CloneTracker[], 1, 0, structure)
module2 = SomaticEvolution.TreeModule(Union{BinaryNode{SimpleTreeCell}, Nothing}[node26, node27, node28], 5.0, [2.0, 4.0], SomaticEvolution.CloneTracker[], 2, 1, structure)
module3 = SomaticEvolution.TreeModule(Union{BinaryNode{SimpleTreeCell}, Nothing}[node29, node30, node23], 5.0, [4.0], SomaticEvolution.CloneTracker[], 3, 2, structure)
population = [module1, module2, module3]
#endregion

@testset "pairwise differences" begin
    @test SomaticEvolution.findMRCA(module1) == mrca1
    @test SomaticEvolution.findMRCA(module2) == mrca2
    @test SomaticEvolution.findMRCA(module3) == mrca3
    @test pairwise_fixed_differences_clonal(population) == (Dict(68 => 1, 102 => 1, 54 => 1), Dict(42 => 1, 32 => 1, 66 => 1))
end

@testset "allelefreq" begin
    @testset "cell subset size" begin
        root = getsingleroot(module1.cells)
        @test SomaticEvolution.cell_subset_size(root, mapreduce(x->x.cells, vcat, population)) == 9
        @test SomaticEvolution.cell_subset_size(root, module1.cells) == 3
        @test SomaticEvolution.cell_subset_size(root, module2.cells) == 3
        @test SomaticEvolution.cell_subset_size(root, module3.cells) == 3
        @test SomaticEvolution.cell_subset_size(root.left.left, module1.cells) == 3
        @test SomaticEvolution.cell_subset_size(root.left.left, module2.cells) == 0
    end
    allcellnodes = mapreduce(x -> x.cells, vcat, population)
    @test SomaticEvolution.cell_subset_size(getsingleroot(allcellnodes), allcellnodes) == length(allcellnodes)
    @test countmap(Vector{Int64}(getallelefreq(module1, 2).*6)) == Dict(1 => 68, 2 => 20, 3 => 42)
    @test countmap(Vector{Int64}(getallelefreq(population, 2).*18)) == Dict(1 => 244, 2 => 106, 3 => 93, 6 => 19, 9 => 3)
end

@testset "module splitting with replacement" begin
    rng = MersenneTwister(12)
    parentmodule = deepcopy(module1)
    parentmodule, newmodule, nextID = SomaticEvolution.sample_new_module_with_replacement!(parentmodule, 2, 1, 
        300, 26, 2, :fixed, rng)
    @test length(parentmodule) == 3
    @test length(newmodule) == 1
    @test nextID == 26+2
end

@testset "simulate module moran" begin

    #check runsimulation function
    rng = MersenneTwister(12)
    input = MultilevelBranchingMoranInput(
            modulesize=4, 
            fixedmu=true, 
            birthrate=0.1, 
            deathrate=0.01,
            moranrate=0.01, 
            clonalmutations=0, 
            tmax=365*4, 
            maxmodules=5,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = runsimulation(SimpleTreeCell, input, rng)
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    @test length(population) == 5

    population = runsimulation(TreeCell, input, rng)
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    @test length(population) == 5
    # @test age(population) < 365 * 4

end





