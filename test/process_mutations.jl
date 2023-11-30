homeostatic_modules = 
    [SomaticEvolution.CellModule(
        Union{Cell, Nothing}[Cell([2, 3, 9, 15], 1, 0), Cell([2, 4], 1, 0), Cell([2, 3, 10], 1, 0), Cell([2, 3, 9, 16], 1, 0)], 
        227.24560639149834,
        [0.0, 166.45694063272188, 218.77527717255302],
        1, 
        0,
        WellMixed()
    )]
growing_modules = [
    SomaticEvolution.CellModule(
        Union{Cell, Nothing}[Cell([1, 6, 8, 11, 13], 1, 0), Cell([1, 6, 8, 12], 1, 0), Cell([1, 6, 8, 11, 14], 1, 0)], 
        216.40410240545967,
        [166.45694063272188], 
        2, 
        1,
        WellMixed()

    )
    SomaticEvolution.CellModule(
        Union{Cell, Nothing}[Cell(Int64[1, 5], 1, 0)], 
        218.77527717255302,
        [218.77527717255302], 
        3, 
        1,
        WellMixed()
    )
]

population = Population(homeostatic_modules, growing_modules, 1.0, 0.0, 1.0, 0.0)

μ = 2
clonalmutations = 3
rng = MersenneTwister(12)

@testset "multilevel" begin
    #gets list of all mutations in the population
    @test Set(SomaticEvolution.get_mutationlist(population)) == 
        Set([1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16])
    @test Set(SomaticEvolution.get_mutationlist(population.homeostatic_modules[1])) ==
        Set([2, 3, 9, 15, 4, 10, 16])
    mutationlist = SomaticEvolution.get_mutationlist(population)

    #check mean is approx correct for poisson distributed mutations
    expandedmutationids = SomaticEvolution.get_expandedmutationids(3, collect(1:100000), clonalmutations, rng, mutationdist=:poisson)
    @test  mean(reduce(vcat,map(length,values(expandedmutationids)))) ≈ 3 atol=0.01
    #check mean is approx correct for geometric distributed mutations
    expandedmutationids = SomaticEvolution.get_expandedmutationids(3, collect(1:100000), clonalmutations, rng, mutationdist=:geometric)
    @test  mean(reduce(vcat,map(length,values(expandedmutationids)))) ≈ 3 atol=0.01

    module1 = deepcopy(population.growing_modules[1])
    expandedmutationids = Dict(1=>[4], 2=>[5,6], 3=>[7], 6=>[8,9,10], 8=>[11], 11=>[], 12=>[], 13=>[12,13,14], 14=>[15,16])
    SomaticEvolution.expandmutations!(module1, expandedmutationids, clonalmutations)
    @test sort(SomaticEvolution.clonal_mutation_ids(module1)) == [1, 2, 3, 4, 8, 9, 10, 11]
    @test sort(module1.cells[1].mutations) == [1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14]
    @test sort(module1.cells[2].mutations) == [1, 2, 3, 4, 8, 9, 10, 11]
    @test sort(module1.cells[3].mutations) == [1, 2, 3, 4, 8, 9, 10, 11, 15, 16]

    expandedmutationids = SomaticEvolution.get_expandedmutationids(μ, mutationlist, clonalmutations, rng, mutationdist=:fixed)
    @test Set(mutationlist) == keys(expandedmutationids)
    @test all(reduce(vcat,map(length,values(expandedmutationids))) .== 2)
    SomaticEvolution.processresults!(population, μ, clonalmutations, rng, mutationdist=:fixed)
    @test clonal_mutations(population.homeostatic_modules[1]) == 2 * 1 + 3
    @test length(population.homeostatic_modules[1].cells[1].mutations) == 2 * 1 + 3 + 2 * 3
    @test length(population.homeostatic_modules[1].cells[2].mutations) == 2 * 1 + 3 + 2 * 1
    @test length(population.growing_modules[2].cells[1].mutations) == 2 * 2 + 3

end