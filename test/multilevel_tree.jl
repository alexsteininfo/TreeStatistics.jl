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
    mutationdist = [:fixed]
    μ = [1]
    modulebranching=:split

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
    population = SomaticEvolution.initialize_population(
        TreeCell, 
        WellMixed, 
        input; 
        rng)

    nextID = SomaticEvolution.getnextID(population)
    nextmoduleID = 2
    @test nextID == 2
    t = 0
    rng = MersenneTwister(1)
    rates = (birthrate=input.birthrate, deathrate=input.deathrate, moranrate=input.moranrate, asymmetricrate=input.asymmetricrate)
    @test SomaticEvolution.getwildtyperates(population) == rates
    transitionrates = SomaticEvolution.get_neutral_transitionrates(
        population, 
        input.branchrate,
        4
    ) 
    #check first update is correct
    #only tranisition with non-zero rate is birth which has rate birthrate=0.1
    @test transitionrates == [0.0, 0.0, 0.1, 0.0, 0.0] 
    population, t, nextID = 
        SomaticEvolution.update_population_neutral!(
            population, 
            transitionrates,
            branchrate, 
            modulesize, 
            branchinitsize,
            modulebranching,
            SomaticEvolution.NoQuiescence(),
            nothing,
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
            mutationdist=:fixed, 
            birthrate=0.1, 
            deathrate=0.01,
            moranrate=0.01, 
            asymmetricrate=0.01, 
            clonalmutations=0, 
            tmax=20*365, 
            maxmodules=10,
            branchrate=3/365, 
            branchinitsize=1, 
            μ=1
    )
    population = runsimulation(SimpleTreeCell, input, rng)
    @test length(population) == 10
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    population = runsimulation(TreeCell, input, rng)
    @test length(population) == 10
    input = MultilevelBranchingInput(
            modulesize=4, 
            mutationdist=:fixed, 
            birthrate=1, 
            deathrate=0,
            moranrate=0.1, 
            clonalmutations=0, 
            tmax=365*100, 
            maxmodules=10,
            branchrate=3/365, 
            branchinitsize=1, 
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
module1 = SomaticEvolution.TreeModule(Union{BinaryNode{SimpleTreeCell}, Nothing}[node19, node24, node25], 5.0, [0.0, 2.0], 1, 0, structure)
module2 = SomaticEvolution.TreeModule(Union{BinaryNode{SimpleTreeCell}, Nothing}[node26, node27, node28], 5.0, [2.0, 4.0], 2, 1, structure)
module3 = SomaticEvolution.TreeModule(Union{BinaryNode{SimpleTreeCell}, Nothing}[node29, node30, node23], 5.0, [4.0], 3, 2, structure)
population = Population([module1, module2, module3], typeof(module1)[], Subclone[Subclone()])
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

@testset "mutation statistics" begin
    @test clonal_mutations.(population) == [42, 32, 66]
    @test all(average_mutations.(population) .≈ [78.0, 106.0, 108.0])
end

@testset "module splitting with replacement" begin
    rng = MersenneTwister(12)
    parentmodule = deepcopy(module1)
    subclones = Subclone[Subclone()]
    parentmodule, newmodule, nextID = 
        SomaticEvolution.sample_new_module_with_replacement!(parentmodule, subclones,  2, 1, 
            300, 26, [2], [:fixed], rng; timedepmutationsonly=false)
    @test length(parentmodule) == 3
    @test length(newmodule) == 1
    @test nextID == 26+2
end

@testset "module splitting without replacement" begin
    #time-dep mutations only
    rng = MersenneTwister(12)
    parentmodule = deepcopy(module1)
    expectedmutations = sort!([round(Int64, 6 * (300 - cell.data.birthtime) + cell.data.mutations) for cell in parentmodule.cells])
    subclones = Subclone[Subclone()]
    parentmodule, newmodule, nextID = 
        SomaticEvolution.sample_new_module_without_replacement!(parentmodule, subclones,  2, 3, 
            300, 26, [2, 6], [:fixed, :fixedtimedep], rng; timedepmutationsonly=true)
    @test length(parentmodule) == 3
    @test length(newmodule) == 3
    @test nextID == 26+6
    @test sort!([cell.parent.data.mutations for cell in parentmodule.cells]) == expectedmutations
    @test sort!([cell.parent.data.mutations for cell in newmodule.cells]) == expectedmutations
    @test [cell.data.mutations for cell in parentmodule.cells] == [0,0,0]
    @test [cell.data.mutations for cell in newmodule.cells] == [0,0,0]

    #all mutations
    rng = MersenneTwister(12)
    parentmodule = deepcopy(module1)
    expectedmutations = sort!([round(Int64, 6 * (300 - cell.data.birthtime) + cell.data.mutations) for cell in parentmodule.cells])
    subclones = Subclone[Subclone()]
    parentmodule, newmodule, nextID = 
        SomaticEvolution.sample_new_module_without_replacement!(parentmodule, subclones,  2, 3, 
            300, 26, [2, 6], [:fixed, :fixedtimedep], rng; timedepmutationsonly=false)
    @test length(parentmodule) == 3
    @test length(newmodule) == 3
    @test nextID == 26+6
    @test sort!([cell.parent.data.mutations for cell in parentmodule.cells]) == expectedmutations
    @test sort!([cell.parent.data.mutations for cell in newmodule.cells]) == expectedmutations
    @test [cell.data.mutations for cell in parentmodule.cells] == [2,2,2]
    @test [cell.data.mutations for cell in newmodule.cells] == [2,2,2]
end

@testset "simulate module moran" begin

    #check runsimulation function
    rng = MersenneTwister(12)
    input = MultilevelBranchingMoranInput(
            modulesize=4, 
            mutationdist=:fixed, 
            birthrate=0.1, 
            deathrate=0.01,
            moranrate=0.01, 
            clonalmutations=0, 
            tmax=365*4, 
            maxmodules=5,
            branchrate=3/365, 
            branchinitsize=1, 
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


@testset "simulate seasonal quiescence" begin

    #check runsimulation function
    rng = MersenneTwister(12)
    quiescence=SeasonalQuiescence(factor=0.5, duration=0.5)
    input = MultilevelBranchingMoranInput(;
            modulesize=4, 
            mutationdist=:fixed, 
            birthrate=100, 
            deathrate=0.0,
            moranrate=100, 
            clonalmutations=0, 
            tmax=4.3, 
            maxmodules=5,
            branchrate=5, 
            branchinitsize=1, 
            μ=1,
            modulebranching=:withoutreplacement_nomutations,
            quiescence
    )
    population = runsimulation(SimpleTreeCell, input, rng)
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    @test length(population) == 5
    # quiescence=SeasonalQuiescence(factor=0.5, duration=0.5)
    # seasonalstate = SomaticEvolution.initializeseason(quiescence)
    # population = Population([mt1], [mt2], 100, 0.0, 100, 20)
    # branchrate = 5
    # @test !seasonalstate.winter
    # @test seasonalstate.nextswitchtime == 0.5
    # transitionrates = SomaticEvolution.get_neutral_transitionrates(population, branchrate, 4, quiescence, seasonalstate)
    # @test transitionrates == Float64[400, 80, 300, 0, 5]
    # t = 0.55
    # seasonswitch, seasonalstate = SomaticEvolution.switchseasons(t, quiescence, seasonalstate)
    # @test seasonswitch
    # @test seasonalstate.winter
    # @test seasonalstate.nextswitchtime == 1.0
    # transitionrates = SomaticEvolution.update_neutral_transitionrates!(transitionrates, population, branchrate, 4, quiescence, seasonalstate)
    # @test transitionrates == Float64[200, 40, 150, 0, 2.5]
end


@testset "simulate stochastic quiescence" begin

    #check runsimulation function
    rng = MersenneTwister(12)
    quiescence=StochasticQuiescence(factor=0.4, onrate=2, offrate=4)
    input = MultilevelBranchingMoranInput(;
            modulesize=4, 
            mutationdist=:fixed, 
            birthrate=100, 
            deathrate=0.0,
            moranrate=100, 
            clonalmutations=0, 
            tmax=4.3, 
            maxmodules=5,
            branchrate=5, 
            branchinitsize=1, 
            μ=1,
            modulebranching=:withoutreplacement_nomutations,
            quiescence
        )
    population = runsimulation(SimpleTreeCell, input, rng).output
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    @test length(population) == 5

    # population = PopulationWithQuiescence([mt1], CellModule{WellMixed}[], [mt2], 100, 0.0, 100, 0.0)
    # transitionrates = SomaticEvolution.get_neutral_transitionrates(population, 5, 4, quiescence)
    # @test transitionrates == Float64[400, 0, 0, 0, 300, 0, 5, 0, 2, 0]
    # SomaticEvolution.quiescenceonupdate!(population, rng) 
    # @test length(population.homeostatic_modules) == 0
    # @test length(population.quiescent_modules) == 1
    # SomaticEvolution.update_neutral_transitionrates!(transitionrates, population, 5, 4, quiescence)
    # @test transitionrates == Float64[0, 160, 0, 0, 300, 0, 0, 2, 0, 4]
end



