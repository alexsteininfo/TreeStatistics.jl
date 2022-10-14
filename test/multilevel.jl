@testset "simulate to tmax" begin

    rng = MersenneTwister(12)
    input = MultilevelBranchingInput(
        modulesize=4, 
        fixedmu=true, 
        b=0.1, 
        d=0.1,
        bdrate=0.1, 
        clonalmutations=0, 
        tmax=400, 
        maxmodules=1000,
        branchrate=3/365, 
        branchfraction=0.2, 
        μ=1
    )
    population = runsimulation(input, rng, :fixedtime)
    @test age(population) < 400
    
    b, d = 0.1, 0.0
    modulesize = 6
    branchinitsize = 2
    branchrate = 1/5
    tmax = 100

    cellmodule = SomaticEvolution.initialize_population(
        modulesize, 
        clonalmutations=0
    )[1]


    rng = MersenneTwister(12)

    cellmodule, newcellmodule =
        SomaticEvolution.module_simulate_to_branching!(
            cellmodule, 
            tmax, 
            b, 
            d, 
            b, 
            branchrate, 
            modulesize, 
            branchinitsize, 
            2, 
            rng
    )
    @test cellmodule.tvec[end] <= tmax
    @test cellmodule.Nvec[end] == modulesize - branchinitsize
    @test newcellmodule.Nvec[end] == branchinitsize
    @test newcellmodule.tvec[end] == cellmodule.tvec[end]

    newtmax = cellmodule.tvec[end] + 0.1
    cellmodule, newcellmodule =
        SomaticEvolution.module_simulate_to_branching!(
            cellmodule, 
            newtmax, 
            b, 
            d, 
            b, 
            branchrate, 
            modulesize, 
            branchinitsize, 
            2, 
            rng
    )
    @test newcellmodule === nothing
end

input = MultilevelBranchingInput(
    modulesize=4, 
    fixedmu=true, 
    b=0.1, 
    d=0,
    bdrate=0.01, 
    clonalmutations=0, 
    tmax=1*365, 
    branchrate=3/365, 
    branchfraction=0.2, 
    μ=1
)

mt1 = SomaticEvolution.CellModule(
    [1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4], 
    [0.0, 14.437821887595856, 38.636412371018324, 38.80067279227201, 39.07905423063021, 56.24183190197407, 61.57555456465887, 93.99343227770015, 174.7383395589789, 177.87796865140686, 187.1205185574646, 190.3156218593258, 222.64188808574013, 235.8640994797343, 242.08245554076214], 
    Cell[Cell([3, 4, 20], 1, 0), Cell([3, 5, 6, 8, 10, 17, 18, 19], 1, 0), Cell([3, 4, 21, 22, 23], 1, 0), Cell([1, 24], 1, 0)], 
    SomaticEvolution.CloneTracker[], 1, 0)
mt2 = SomaticEvolution.CellModule(
    [1, 2, 3, 4, 4, 4, 3], 
    [187.1205185574646, 207.5419622150601, 209.92705453342273, 216.2145518957684, 243.24240821372297, 252.31115142542512, 257.22794191422554], 
    Cell[Cell([1, 13, 25, 26, 27], 1, 0), Cell([1, 13], 1, 0), Cell([1, 13, 25, 28, 29], 1, 0)], 
    SomaticEvolution.CloneTracker[], 2, 1)
mt3 = SomaticEvolution.CellModule(
    [1], 
    [257.22794191422554], 
    Cell[Cell([1], 1, 0)], SomaticEvolution.CloneTracker[], 3, 2)

@testset "module sampling" begin
    rng = MersenneTwister(12)
    cellmodule = deepcopy(mt1)
    cellmodule, newcellmodule = 
        SomaticEvolution.sample_new_module!(cellmodule, 2, 1, 1.0, rng)
    @test length(cellmodule) == 3
    @test length(newcellmodule) == 1
    @test newcellmodule.id == 2
    @test newcellmodule.Nvec == [1]
    @test newcellmodule.tvec == [1.0]
    cellmodule = deepcopy(mt1)
    nextmoduleID = 2
    population, nextmoduleID = SomaticEvolution.modulebranchingupdate!([cellmodule], nextmoduleID, 4, 1, 2.0, rng)
    @test nextmoduleID == 3
    @test moduleid.(population) == [1,2]
end
    

population = SomaticEvolution.MultiSimulation(input, [mt1, mt2, mt3])

#check that statistics are calculated correctly
@testset "mutation statistics" begin
    @test mutations_per_cell(mt1) == [3, 8, 5, 2]
    @test average_mutations_per_module(population) ≈ [4.5, 4, 1]
    @test all(average_mutations(population, true) .≈ (3.875, 5.267857142857143))
    @test clonal_mutations(population) == [0, 2, 1]
    @test SomaticEvolution.clonal_mutation_ids(population) == [[], [1, 13], [1]]
    @test pairwise_fixed_differences_matrix(population, diagonals=true) == [0 0 0; 2 2 0; 1 1 1]
    @test pairwise_fixed_differences_matrix(population, [2,3]) == [0 0; 1 0]
    @test pairwise_fixed_differences_clonal(population) == (Dict{Int64, Int64}(1=>2, 2=>1), Dict{Int64, Int64}(0=>1, 2=>1, 1=>1))
    @test pairwise_fixed_differences_clonal(population, [2,3]) == (Dict{Int64, Int64}(1=>1), Dict{Int64, Int64}(2=>1, 1=>1))
    @test shared_fixed_mutations(population) == Dict{Int64, Int64}(1=>1, 2=>1)
    @test shared_fixed_mutations(population, [2,3]) == Dict{Int64, Int64}(1=>1, 2=>1)
    @test all(pairwise_fixed_differences_statistics(population) .≈ (1.3333333333333333,0.33333333333333333,1,1))
    @test newmoduletimes(population) ≈ [0.0, 187.1205185574646, 257.22794191422554]
    @test all(cellpopulationsize(population, 50) .≈ ([0,50,100,150,200,250],[1, 4, 4, 4, 5, 8]))
    @test all(meanmodulesize(population, 50) .≈ ([0,50,100,150,200,250], [1.0, 4.0, 4.0, 4.0, 2.5, 4.0]))
    @test all(numbermodules(population, 50) .≈ ([0,50,100,150,200,250], [1, 1, 1, 1, 2, 2]))
end

@testset "updates" begin
    rng = MersenneTwister(12)
    population = [deepcopy(mt1), deepcopy(mt1), deepcopy(mt3)]

    #kills only cell in module 3 => population size goes to 2 
    SomaticEvolution.deathupdate!(population, 4, 3, 1, :fixed, rng)
    @test length(population) == 2

    #one module splits into two modules of length two
    SomaticEvolution.modulebranchingupdate!(population, 2, 4, 2, 1.0, rng)
    @test length(population) == 3
    @test sum(map(x -> length(x.cells), population)) == 8
    @test sum(map(x -> x.Nvec[end], population)) == 8

    #birth update one module gets an extra cell
    SomaticEvolution.birthupdate!(population, 4, 1.0, 2, 1, :fixed, rng)
    @test sum(map(x -> length(x.cells), population)) == 9
    @test sum(map(x -> x.Nvec[end], population)) == 9

    #moran update population size stays the same
    SomaticEvolution.moranupdate!(population, 4, 1.0, 3, 1, :fixed, rng)
    @test sum(map(x -> length(x.cells), population)) == 9
    @test sum(map(x -> x.Nvec[end], population)) == 9

    N0 = population[1].Nvec[end]
    @test N0 == length(population[1].cells)
    L = length(population[1].Nvec)
    @test L == length(population[1].tvec)
    SomaticEvolution.updatemodulehistory!(population[1], 1, 20.0)
    @test population[1].tvec[end] == 20.0
    @test population[1].Nvec[end] == 1 + N0
    @test length(population[1].Nvec) == length(population[1].tvec) == L + 1
    @test_throws ErrorException SomaticEvolution.check_module_number(population, 1)


end

@testset "simulate to fixed size" begin
    b = 0.1
    d = 0.0
    bdrate = 0.01
    modulesize = 4
    branchinitsize = 1
    branchrate = 3/365
    tmax = 1000
    maxmodules = 5
    population = SomaticEvolution.initialize_population(modulesize, clonalmutations=0)
    mutID = SomaticEvolution.getnextID(population[1].cells) #get id of next mutation
    t = 0
    rng = MersenneTwister(1)

    #check first update is correct
    transitionrates = SomaticEvolution.get_transitionrates(
        population, 
        b, 
        d, 
        bdrate, 
        branchrate, 
        modulesize
    ) 
    #only tranisition with non-zero rate is birth which has rate b=0.1
    @test transitionrates == [0.0, 0.1, 0.0, 0.0] 
    population, t, mutID, nextmoduleID = 
        SomaticEvolution.update_population!(
            population, 
            b, 
            d, 
            bdrate, 
            branchrate, 
            modulesize, 
            branchinitsize,
            t,
            mutID, 
            2,
            1,
            :fixed, 
            tmax,
            maxmodules,
            rng
    )
    #popsize is now 2
    @test population[1].Nvec[end] == length(population[1].cells) == 2
    @test length(population[1].Nvec) === length(population[1].tvec) == 2
    #each cell has one mutation
    @test length(population[1].cells[1].mutations) == length(population[1].cells[2].mutations) == 1

    #check transition rates for population of modules m1, m2, m3
    transitionrates = SomaticEvolution.get_transitionrates(
        [mt1, mt2, mt3], 
        b, 
        d, 
        bdrate, 
        branchrate, 
        modulesize
    ) 
    @test all(transitionrates .≈ [0.04, 0.4, 0.0, 3/365])

    #check transitionrates for non-zero death rate
    transitionrates = SomaticEvolution.get_transitionrates(
        [mt1, mt2, mt3], 
        b, 
        0.01, 
        bdrate, 
        branchrate, 
        modulesize
    ) 
    @test all(transitionrates .≈ [0.04, 0.4, 0.04, 3/365])

    #check only homeostatic modules are chosen for Moran and module branching and that prob of choosing cells is equal
    @test all(SomaticEvolution.choose_homeostaticmodule([mt1,mt2,mt3], 4, rng)==mt1 for i in 1:10)
    @test isapprox(mean(reduce(hcat, [[v for v in SomaticEvolution.choose_homeostaticmodule_cells([mt1,mt2,mt3], 4, rng)[2:3]] for i in 1:10000]), dims=2), [2.5, 2.5], atol=0.1)
    #check prob of choosing non-homeostatic modules is proportional to number of cells
    @test isapprox(mean(length(SomaticEvolution.choose_growingmodule_cell([mt1,mt2,mt3], 4, rng)[1].cells) for i in 1:100000), 2.5, atol=0.1)

    #check runsimulation function
    rng = MersenneTwister(12)
    input = MultilevelBranchingInput(
            modulesize=4, 
            fixedmu=true, 
            b=0.1, 
            d=0.01,
            bdrate=0.01, 
            clonalmutations=0, 
            tmax=20*365, 
            maxmodules=10,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = runsimulation(input, rng, :normal)
    @test length(population) == 10
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    input = MultilevelBranchingInput(
            modulesize=4, 
            fixedmu=true, 
            b=0.1, 
            d=0.1,
            bdrate=0.1, 
            clonalmutations=0, 
            tmax=365, 
            maxmodules=1000,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = runsimulation(input, rng, :normal)
    @test age(population) <= 365

end

mt1 = SomaticEvolution.CellModule(
    [1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4], 
    [0.0, 14.437821887595856, 38.636412371018324, 38.80067279227201, 39.07905423063021, 56.24183190197407, 61.57555456465887, 93.99343227770015, 174.7383395589789, 177.87796865140686, 187.1205185574646, 190.3156218593258, 222.64188808574013, 235.8640994797343, 242.08245554076214], 
    Cell[Cell([3, 4, 20], 1, 0), Cell([3, 5, 6, 8, 10, 17, 18, 19], 1, 0), Cell([3, 4, 21, 22, 23], 1, 0), Cell([1, 24], 1, 0)], 
    SomaticEvolution.CloneTracker[], 1, 0)
mt2 = SomaticEvolution.CellModule(
    [1, 2, 3, 4, 4, 4, 3], 
    [187.1205185574646, 207.5419622150601, 209.92705453342273, 216.2145518957684, 243.24240821372297, 252.31115142542512, 257.22794191422554], 
    Cell[Cell([1, 13, 25, 26, 27], 1, 0), Cell([1, 13], 1, 0), Cell([1, 13, 25, 28, 29], 1, 0)], 
    SomaticEvolution.CloneTracker[], 2, 1)
mt3 = SomaticEvolution.CellModule(
    [1], 
    [257.22794191422554], 
    Cell[Cell([1], 1, 0)], SomaticEvolution.CloneTracker[], 3, 2)

@testset "moran updates" begin
    rng = MersenneTwister(12)
    population = [deepcopy(mt1), deepcopy(mt2), deepcopy(mt3)]
    transitionrates = SomaticEvolution.get_transitionrates(
        population, 
        1.0, 
        0, 
        1.0, 
        0.1, 
        4
    ) 
    @test transitionrates == [4.0, 4.0, 0.0, 0.1]
    SomaticEvolution.modulemoranupdate!(population, 4, 4, 2, 5, 1, :fixed, rng)
    length(population) == 3


end


@testset "simulate module moran" begin

    #check runsimulation function
    rng = MersenneTwister(12)
    input = MultilevelBranchingMoranInput(
            modulesize=4, 
            fixedmu=true, 
            b=0.1, 
            d=0.01,
            bdrate=0.01, 
            clonalmutations=0, 
            tmax=365*4, 
            maxmodules=5,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = runsimulation(input, rng)
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    @test length(population) == 5
    # @test age(population) < 365 * 4

end
