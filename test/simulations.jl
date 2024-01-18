@testset "neutral runsimulation" begin
    rng = MersenneTwister(100)
    input = BranchingInput(
        Nmax=10, 
        mutationdist=:poisson, 
        birthrate=1, 
        deathrate=0.0,
        clonalmutations=0, 
        μ=1
    )
    simulation = runsimulation(Cell, input, rng)
    @test length(simulation.output) == 10

    tmax=10
    input = MoranInput(
        N=10, 
        tmax=tmax,
        mutationdist=:poisson, 
        moranrate=1.0,
        clonalmutations=0, 
        μ=1
    )
    simulation = runsimulation(Cell, input, rng)
    @test all(length(simulation.output) .== 10)
    @test age(simulation) == simulation.output.singlemodule.t
    @test age(simulation) <=tmax
    tmax=10
    input = BranchingMoranInput(
        Nmax=10, 
        tmax=tmax,
        mutationdist=:poisson, 
        birthrate=1, 
        deathrate=0.0,
        clonalmutations=0, 
        μ=1
    )
    simulation = runsimulation(Cell, input, rng)
    @test length(simulation.output) == 10
    @test age(simulation) <= tmax
end

mt1 = SomaticEvolution.CellModule(
    Union{Cell, Nothing}[Cell([1, 2], 1, 0.3), Cell([1, 3], 1, 0.4), Cell([1, 4], 1, 0.4)], 
    0.4,
    [0.0],
    1, 0, WellMixed())

@testset "updates" begin
    rng = MersenneTwister(12)
    population = SinglelevelPopulation(
        deepcopy(mt1), Subclone[Subclone(1, 0, 0.0, 3, 1.0, 0.1, 1.0, 0.1)]
    )
    SomaticEvolution.celldivision!(population.singlemodule, population.subclones, 1, 0.5, 5, [1], [:fixed], rng)  
    @test length(population) == 4
    @test getsubclonesizes(population) == [4]
    SomaticEvolution.cellmutation!(population.singlemodule, population.subclones, 0.5, population.singlemodule.cells[1], 0.5)
    @test getclonetype(population.singlemodule.cells[1]) == 2
    @test length(population.subclones) == 2
    @test population.subclones[2].birthrate ≈ 1.5
    @test population.subclones[2].deathrate ≈ 0.1
    @test population.subclones[2].moranrate ≈ 1.5
    @test population.subclones[2].asymmetricrate ≈ 0.15
    @test getsubclonesizes(population) == [3,1]
    SomaticEvolution.moranupdate!(
        population, 
        SelectionPredefined(Float64[0.5], Float64[0.5]),
        SomaticEvolution.getmoranrates(population.subclones),
        maximum(SomaticEvolution.getmoranrates(population.subclones)),
        4, 7, 2, 2, 0.51, [1], [:fixed], false, rng
    )
    @test length(population) == 4
end

@testset "rates" begin
#check birth and death rates are calculated correctly
    subclone = SomaticEvolution.Subclone(1, 0, 0.0, 1, 2.0, 1.0, 1.0, 2.0)
    @test SomaticEvolution.getwildtyperates([subclone]) == (birthrate=2.0, deathrate=1.0, moranrate=1.0, asymmetricrate=2.0)
    newrates = SomaticEvolution.get_newsubclone_rates(SomaticEvolution.getwildtyperates([subclone]), 0.5)
    @test newrates == (birthrate=3.0, deathrate=1.0, moranrate=1.5, asymmetricrate=3.0)
end

@testset "selection runsimulation" begin
    rng = MersenneTwister(100)
    input = BranchingInput(
        Nmax=10, 
        mutationdist=:poisson, 
        birthrate=1, 
        deathrate=0.0,
        clonalmutations=0, 
        μ=1
    )
    selection = SelectionPredefined([0.5],[3])
    simulation = runsimulation(Cell, input, selection, rng)
    @test length(simulation.output) == 10
    @test sum(getsubclonesizes(simulation)) == 10
    @test getsubclonesizes(simulation) == counts(SomaticEvolution.getclonetype.(simulation.output.singlemodule.cells))


    tmax=10
    input = MoranInput(
        N=10, 
        tmax=tmax,
        mutationdist=:poisson, 
        moranrate=1.0,
        clonalmutations=0, 
        μ=1
    )
    simulation = runsimulation(Cell, input, rng)
    @test all(length(simulation.output) .== 10)
    @test age(simulation) == simulation.output.singlemodule.t
    @test age(simulation) <= tmax
    tmax=10
    input = BranchingMoranInput(
        Nmax=10, 
        tmax=tmax,
        mutationdist=:poisson, 
        birthrate=10, 
        deathrate=0.0,
        clonalmutations=0, 
        μ=1
    )
    selection = SelectionPredefined([0.5, 1.0], [0.1, 7.0])
    simulation = runsimulation(Cell, input, selection, rng)
    @test length(simulation.output) == 10
    @test sum(getsubclonesizes(simulation)) == 10
    subclone_by_cell = SomaticEvolution.getclonetype.(simulation.output.singlemodule.cells)
    @test getsubclonesizes(simulation) == counts(subclone_by_cell, 1:length(simulation.output.subclones))
    @test age(simulation) <= tmax
    @test simulation.output.subclones[2].mutationtime ≈ 0.1 atol=0.1
    @test simulation.output.subclones[3].mutationtime ≈ 7.0 atol=1
end 


@testset "mutations" begin
    input = BranchingInput(
        clonalmutations=0, 
        Nmax=1, 
        birthrate=1, 
        deathrate=0,
        μ=[1, 3],
        mutationdist=[:fixed, :fixedtimedep]
    )
    population = SomaticEvolution.initialize_population(Cell, WellMixed, input; rng)
    treemodule = population.singlemodule
    subclones = population.subclones
    nextID = 1
    @test treemodule.cells[1].mutations == Int64[]
    treemodule, subclones, nextID = SomaticEvolution.celldivision!(treemodule, subclones, 1, 1.0, nextID, input.μ, input.mutationdist, rng)
    @test treemodule.cells[1].mutations == Int64[1,3,4,5]
    @test treemodule.cells[2].mutations == Int64[2,3,4,5]
    treemodule, subclones, nextID = SomaticEvolution.celldivision!(treemodule, subclones, 1, 3.0, nextID, input.μ, input.mutationdist, rng)
    @test treemodule.cells[1].mutations == Int64[1,3,4,5,6,8,9,10,11,12,13]
    @test treemodule.cells[2].mutations == Int64[2,3,4,5]
    @test treemodule.cells[3].mutations == Int64[1,3,4,5,7,8,9,10,11,12,13]
    SomaticEvolution.final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=4.0)
    @test treemodule.cells[1].mutations == Int64[1,3,4,5,6,8,9,10,11,12,13,14,15,16]
    @test treemodule.cells[2].mutations == Int64[2,3,4,5,17,18,19,20,21,22,23,24,25]
    @test treemodule.cells[3].mutations == Int64[1,3,4,5,7,8,9,10,11,12,13,26,27,28]
end
