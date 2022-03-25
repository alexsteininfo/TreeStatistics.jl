populationtracker = [
    SomaticEvolution.ModuleTracker(
        [1, 2, 3, 4, 4, 3, 4, 3, 4], 
        [0.0, 16.892494137270344, 68.68837372525104, 73.9360343780634, 95.06744895741765, 166.45694063272188, 186.08040859501872, 218.77527717255302, 227.24560639149834], 
        Cell[Cell([2, 3, 9, 15], 1), Cell([2, 4], 1), Cell([2, 3, 10], 1), Cell([2, 3, 9, 16], 1)], 
        SomaticEvolution.CloneTracker[], 
        1, 
        0
    )
    SomaticEvolution.ModuleTracker(
        [1, 2, 3], 
        [166.45694063272188, 215.96281393553738, 216.40410240545967], 
        Cell[Cell([1, 6, 8, 11, 13], 1), Cell([1, 6, 8, 12], 1), Cell([1, 6, 8, 11, 14], 1)], 
        SomaticEvolution.CloneTracker[], 
        2, 
        1
    )
    SomaticEvolution.ModuleTracker(
        [1], 
        [218.77527717255302], 
        Cell[Cell(Int64[1, 5], 1)], 
        SomaticEvolution.CloneTracker[], 
        3, 
        1
    )
]

μ = 2
clonalmutations = 3
rng = MersenneTwister(12)

@testset "multilevel" begin
    #gets list of all mutations in the population
    @test Set(SomaticEvolution.get_mutationlist(populationtracker)) == 
        Set([1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16])
    @test Set(SomaticEvolution.get_mutationlist(populationtracker[1])) ==
        Set([2, 3, 9, 15, 4, 10, 16])
    mutationlist = SomaticEvolution.get_mutationlist(populationtracker)

    #check mean is approx correct for poisson distributed mutations
    expandedmutationids = SomaticEvolution.get_expandedmutationids(3, collect(1:100000), clonalmutations, rng, fixedmu=false)
    @test  mean(reduce(vcat,map(length,values(expandedmutationids)))) ≈ 3 atol=0.01

    module1 = deepcopy(populationtracker[2])
    expandedmutationids = Dict(1=>[4], 2=>[5,6], 3=>[7], 6=>[8,9,10], 8=>[11], 11=>[], 12=>[], 13=>[12,13,14], 14=>[15,16])
    SomaticEvolution.expandmutations!(module1, expandedmutationids, clonalmutations)
    @test sort(clonal_mutation_ids(module1)) == [1, 2, 3, 4, 8, 9, 10, 11]
    @test sort(module1.cells[1].mutations) == [1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14]
    @test sort(module1.cells[2].mutations) == [1, 2, 3, 4, 8, 9, 10, 11]
    @test sort(module1.cells[3].mutations) == [1, 2, 3, 4, 8, 9, 10, 11, 15, 16]

    expandedmutationids = SomaticEvolution.get_expandedmutationids(μ, mutationlist, clonalmutations, rng, fixedmu=true)
    @test Set(mutationlist) == keys(expandedmutationids)
    @test all(reduce(vcat,map(length,values(expandedmutationids))) .== 2)
    SomaticEvolution.processresults!(populationtracker, μ, clonalmutations, rng, fixedmu=true)
    @test clonal_mutations(populationtracker[1]) == 2 * 1 + 3
    @test length(populationtracker[1].cells[1].mutations) == 2 * 1 + 3 + 2 * 3
    @test length(populationtracker[1].cells[2].mutations) == 2 * 1 + 3 + 2 * 1
    @test length(populationtracker[3].cells[1].mutations) == 2 * 2 + 3

end

SomaticEvolution.ModuleTracker(
    [1, 2, 3, 4, 4, 3, 4, 3, 4], 
    [0.0, 16.892494137270344, 68.68837372525104, 73.9360343780634, 95.06744895741765, 166.45694063272188, 186.08040859501872, 218.77527717255302, 227.24560639149834], 
    Cell[Cell([2, 3, 9, 15], 1), Cell([2, 4], 1), Cell([2, 3, 10], 1), Cell([2, 3, 9, 16], 1)], 
    SomaticEvolution.CloneTracker[], 
    1, 
    0
)

@testset "singlesim" begin
    
end