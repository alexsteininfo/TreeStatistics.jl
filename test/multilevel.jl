input = MultilevelBranchingInput(;
    modulesize=4, 
    mutationdist=:fixed, 
    birthrate=0.1, 
    deathrate=0,
    moranrate=0.01, 
    clonalmutations=0, 
    tmax=1*365, 
    branchrate=3/365, 
    branchinitsize=1, 
    μ=1
)

mt1 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([3, 4, 20], 1, 0), Cell([3, 5, 6, 8, 10, 17, 18, 19], 1, 0), Cell([3, 4, 21, 22, 23], 1, 0), Cell([1, 24], 1, 0)], 
    242.1,
    [0.0, 187.12],
    1, 0, WellMixed())
mt2 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([1, 13, 25, 26, 27], 1, 0), Cell([1, 13], 1, 0), Cell([1, 13, 25, 28, 29], 1, 0)], 
    257.2, 
    [257.2],
    2, 1, WellMixed())
mt3 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([1], 1, 0)],
    257.2,
    [257.2], 
    3, 2, WellMixed())

population1 = Population([deepcopy(mt1)], CellModule{WellMixed}[], [SomaticEvolution.Subclone(;size=4)])

@testset "module sampling" begin
    rng = MersenneTwister(12)
    cellmodule = deepcopy(mt1)
    cellmodule, newcellmodule = 
        SomaticEvolution.sample_new_module_split!(cellmodule, 2, 1, 1.0, rng)
    @test length(cellmodule) == 3
    @test length(newcellmodule) == 1
    @test newcellmodule.id == 2
    @test newcellmodule.t == 1.0
    cellmodule = deepcopy(mt1)
    nextmoduleID = 2
    population, nextmoduleID = SomaticEvolution.modulebranchingupdate!(population1, nextmoduleID, 1, 2.0, rng)
    @test nextmoduleID == 3
    @test moduleid.(population1.growing_modules) == [2, 1]
    @test length(population1.homeostatic_modules) == 0
end
 
rng = Random.MersenneTwister(12)
population2 = Population([deepcopy(mt1)], deepcopy.([mt2, mt3]), [SomaticEvolution.Subclone(;size=8)])
simulation = SomaticEvolution.Simulation(input, population2)

#check that statistics are calculated correctly
@testset "mutation statistics" begin
    @test mutations_per_cell(mt1) == [3, 8, 5, 2]
    @test average_mutations_per_module(simulation) ≈ [4.5, 4, 1]
    @test all(average_mutations(simulation, true) .≈ (3.875, 5.267857142857143))
    @test clonal_mutations(simulation) == [0, 2, 1]
    @test SomaticEvolution.clonal_mutation_ids(simulation) == [[], [1, 13], [1]]
    @test pairwise_fixed_differences_matrix(simulation, diagonals=true) == [0 0 0; 2 2 0; 1 1 1]
    @test pairwise_fixed_differences_matrix(simulation, [2,3]) == [0 0; 1 0]
    @test pairwise_fixed_differences_clonal(simulation) == (Dict{Int64, Int64}(1=>2, 2=>1), Dict{Int64, Int64}(0=>1, 2=>1, 1=>1))
    @test pairwise_fixed_differences_clonal(simulation, [2,3]) == (Dict{Int64, Int64}(1=>1), Dict{Int64, Int64}(2=>1, 1=>1))
    @test shared_fixed_mutations(simulation) == Dict{Int64, Int64}(1=>1, 2=>1)
    @test shared_fixed_mutations(simulation, [2,3]) == Dict{Int64, Int64}(1=>1, 2=>1)
    @test all(pairwise_fixed_differences_statistics(simulation) .≈ (1.3333333333333333,0.33333333333333333,1,1))
    @test countmap(Vector{Int64}(getallelefreq(mt1, 2).*8)) == Dict(1 => 13, 2 => 1, 3 => 1)
    @test countmap(Vector{Int64}(getallelefreq(simulation).*16)) == Dict(1 => 16, 2 => 2, 3 => 2, 5 => 1)
    @test countmap(Vector{Int64}(getallelefreq(simulation, 2).*6)) == Dict(1 => 4, 2 => 1, 3 => 2)

end

@testset "updates" begin
    rng = MersenneTwister(12)
    population = Population([deepcopy(mt1)], [deepcopy(mt3)], [SomaticEvolution.Subclone(;size=8)])

    #moran update population size stays the same
    @test sum(map(x -> length(x.cells), population)) == 5
    SomaticEvolution.moranupdate!(population, 4, 1.0, 3, [1], [:fixed], rng)
    @test sum(map(x -> length(x.cells), population)) == 5

    #kills only cell in module mt3 => population size goes to 1
    SomaticEvolution.deathupdate!(population, 1.0, [1], [:fixed], rng)
    @test length(population) == 1

    #one module splits into two modules of length two
    SomaticEvolution.modulebranchingupdate!(population, 4, 2, 1.0, rng)
    @test length(population) == 2
    @test sum(map(x -> length(x.cells), population)) == 4
    @test sum(map(x -> length(x), population)) == 4

    #birth update one module gets an extra cell
    SomaticEvolution.birthupdate!(population, 4, 1.0, 2, [1], [:fixed], rng)
    @test sum(map(x -> length(x.cells), population)) == 5
    @test sum(map(x -> length(x), population)) == 5

end

modulestructure = WellMixed()
mt1 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([3, 4, 20], 1, 0), Cell([3, 5, 6, 8, 10, 17, 18, 19], 1, 0), Cell([3, 4, 21, 22, 23], 1, 0), Cell([1, 24], 1, 0)], 
    242.1,
    [0.0, 187.12],
    1, 0, modulestructure)
mt2 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([1, 13, 25, 26, 27], 1, 0), Cell([1, 13], 1, 0), Cell([1, 13, 25, 28, 29], 1, 0)], 
    257.2, 
    [257.2],
    2, 1, modulestructure)
mt3 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([1], 1, 0)],
    257.2,
    [257.2], 
    3, 2, modulestructure)

@testset "module splitting with replacement" begin
    rng = MersenneTwister(12)
    parentmodule = deepcopy(mt1)
    subclones = Subclone[Subclone()]
    parentmodule, newmodule, nextID = SomaticEvolution.sample_new_module_with_replacement!(parentmodule, subclones, 2, 1, 
        300, 25, [2], [:fixed], rng; timedepmutationsonly=false)
    @test length(parentmodule) == 4
    @test length(newmodule) == 1
    @test nextID == 25+4
end
@testset "module splitting without replacement" begin
    rng = MersenneTwister(12)
    parentmodule = deepcopy(mt1)
    subclones = Subclone[Subclone()]
    parentmodule, newmodule, nextID = 
        SomaticEvolution.sample_new_module_without_replacement!(
            parentmodule, subclones, 2, 2, 300, 25, [2], [:fixed], rng; timedepmutationsonly=false)
    @test length(parentmodule) == 4
    @test length(newmodule) == 2
    @test length(newmodule.cells) == 2 
    @test newmodule.cells[1] != newmodule.cells[2]
    @test nextID == 25+8

    rng = MersenneTwister(12)
    parentmodule = deepcopy(mt1)
    subclones = Subclone[Subclone()]
    expectedtimedepmutations = sum([round(Int64, 0.01 * (300 - cell.birthtime)) for cell in parentmodule.cells])
    parentmodule, newmodule, nextID = 
        SomaticEvolution.sample_new_module_without_replacement!(
            parentmodule, subclones, 2, 4, 300, 25, [2, 0.01], [:fixed, :fixedtimedep], rng; timedepmutationsonly=true)
    @test length(parentmodule) == 4
    @test length(newmodule) == 4
    @test length(newmodule.cells) == 4 
    @test newmodule.cells[1] != newmodule.cells[2]
    @test nextID == 25+expectedtimedepmutations
end

@testset "module splitting without replacement no mutations" begin
    rng = MersenneTwister(12)
    parentmodule = deepcopy(mt1)
    subclones = Subclone[Subclone()]
    parentmodule, newmodule, nextID = 
        SomaticEvolution.newmoduleformation!(parentmodule, subclones, 5, 2, 2, rng; 
            modulebranching=:withoutreplacement_nomutations, nextID=25, μ=[1], mutationdist=[:poisson])
    @test length(parentmodule) == 4
    @test length(newmodule) == 2
    @test newmodule.cells[1] != newmodule.cells[2]
    @test nextID == 25
    in.(map(x->x.mutations, newmodule.cells), (map(x->x.mutations, parentmodule.cells), ))
end

@testset "moran updates" begin
    rng = MersenneTwister(12)
    population = Population(
        [deepcopy(mt1)], 
        [deepcopy(mt2), deepcopy(mt3)], 
        Subclone[Subclone(;birthrate=1.0, deathrate=0.5, moranrate=0.2, asymmetricrate=1.0)]
    )
    @test SomaticEvolution.getwildtyperates(population) == (birthrate=1.0, deathrate=0.5, moranrate=0.2, asymmetricrate=1.0)
    branchrate = 10.0
    modulesize = 4
    transitionrates = SomaticEvolution.get_neutral_transitionrates(
        population, 
        branchrate,
        4
    ) 
    @test transitionrates == [0.8, 4.0, 4.0, 2.0, 10.0]
    
    population, = SomaticEvolution.moranupdate!(population, 4, 1, 30, [1], [:poisson], rng; moranincludeself=false)
    @test length(mt1) == 4

    population, = SomaticEvolution.asymmetricupdate!(population, 4, 2, 32, [1], [:poisson], rng)
    @test length(mt1) == 4

    SomaticEvolution.modulemoranupdate!(population, 4, 2, 5, rng)
    @test length(population) == 3

    #test moran update not including self
    population = Population(
        [SomaticEvolution.CellModule(
            Union{Cell, Nothing}[Cell([1], 1, 0), Cell([2], 1, 0)], 2.0, [0.0], 
            1, 0, WellMixed())],
        SomaticEvolution.CellModule{WellMixed}[],
        Subclone[Subclone()])
    
    nextID = 3
    for i in 1:10
        chosenmoduleid, dividecellid, deadcellid = SomaticEvolution.choose_homeostaticmodule_cells(population, rng; twocells=true, moranincludeself=false)
        @test dividecellid != deadcellid
        population, nextID = SomaticEvolution.moranupdate!(population, 2, 1, nextID, [1], [:fixed], rng; moranincludeself=false)
        @test population[1].cells[1].mutations[1:end-1] == population[1].cells[2].mutations[1:end-1]
    end

    #test moran update not including self with NonSpatialTreeModule{SimpleTreeCell}
    input = MultilevelBranchingMoranInput(
        modulesize=2, 
        mutationdist=:fixed, 
        birthrate=100, 
        deathrate=0,
        moranrate=100, 
        clonalmutations=0, 
        tmax=1, 
        maxmodules=1,
        branchrate=5, 
        branchinitsize=1, 
        μ=1,
        moranincludeself=false
    )
    population = runsimulation(SimpleTreeCell, input, rng).output
    nextID = maximum(cell.data.id for cell in population[1].cells) + 1
    for i in 1:10
        chosenmoduleid, dividecellid, deadcellid = SomaticEvolution.choose_homeostaticmodule_cells(population, rng; twocells=true, moranincludeself=false)
        @test dividecellid != deadcellid
        population, nextID = SomaticEvolution.moranupdate!(population, 2, 1, nextID, [1], [:fixed], rng; moranincludeself=false)
        @test population[1].cells[1].data.birthtime == population[1].cells[2].data.birthtime
    end


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
            μ=1,
            modulebranching=:withoutreplacement_nomutations
    )
    population = runsimulation(input, rng)
    @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
    @test length(population) == 5
    # @test age(population) < 365 * 4

end
