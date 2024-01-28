@testset "population initialization" begin
    rng = Random.MersenneTwister(12)
    clonalmutations = 50
    initmodule = SomaticEvolution.initialize(Cell, WellMixed, clonalmutations, 1)
    @test length(initmodule) == 1
    @test initmodule.t == 0.0
    @test length(initmodule.cells) == 1
    @test length(initmodule.cells[1].mutations) == clonalmutations
    @test initmodule.cells[1].clonetype == 1
    @test initmodule.id == 1
    birthrate=1.0
    deathrate=0.0
    moranrate=2.0
    asymmetricrate=0.5
    input = MultilevelMoranInput(;clonalmutations, modulesize=5, maxmodules=2, birthrate,
        deathrate, moranrate, asymmetricrate)
    population = SomaticEvolution.initialize_population(Cell, WellMixed, input; rng)
    @test length(population.subclones) == 1
    @test SomaticEvolution.getwildtyperates(population) == (;birthrate, deathrate, moranrate, asymmetricrate)
    @test length(population) == 2
    @test length(population.growing_modules) == 2
    @test length(population.homeostatic_modules) == 0
    @test length.(population.growing_modules) == [1,1]
end

@testset "module sampling" begin
        cells = Cell[]
        push!(cells, Cell([1,2,3,4], 1, 0.2, 0.2, 1, 0))
        push!(cells, Cell([1,2,3], 1, 0.2, 0.2, 2, 1))
        id = 1
        parentid = 0
        inittime = 5
        initmodule = SomaticEvolution.new_module_from_cells(cells, inittime, [inittime], id,
            parentid, WellMixed())
        @test initmodule.branchtimes[1] == inittime
        @test length(initmodule.cells) == length(cells)
        @test initmodule.cells[1].mutations == [1,2,3,4]
        @test SomaticEvolution.getclonetype(initmodule.cells[1]) == SomaticEvolution.getclonetype(initmodule.cells[2]) == 1
        @test initmodule.id == id
end
