@testset "simulate to fixed size" begin
    b = 0.1
    d = 0.0
    bdrate = 0.01
    modulesize = 4
    branchinitsize = 1
    branchrate = 0.01
    maxtime = 1000
    maxmodules = 5
    mutationdist = :fixed
    μ = 1

    input = MultilevelBranchingInput(;b, d, bdrate, modulesize, branchinitsize, branchrate, maxtime, maxmodules, μ, mutationdist)

    rng = MersenneTwister(1)
    populationtracker = SomaticEvolution.initialize_population(TreeCell, input, rng)
    nextID = SomaticEvolution.getnextID(populationtracker) #get id of next cell
    nextmoduleID = 2
    @test nextID == 2
    t = 0
    rng = MersenneTwister(1)

    #check first update is correct
    transitionrates = SomaticEvolution.get_transitionrates(
        populationtracker, 
        b, 
        d, 
        bdrate, 
        branchrate, 
        modulesize
    ) 
    #only tranisition with non-zero rate is birth which has rate b=0.1
    @test transitionrates == [0.0, 0.1, 0.0, 0.0] 
    populationtracker, t, nextID = 
        SomaticEvolution.update_population!(
            populationtracker, 
            b, 
            d, 
            bdrate, 
            branchrate, 
            modulesize, 
            branchinitsize,
            t,
            nextID, 
            nextmoduleID,
            μ,
            mutationdist, 
            maxtime,
            maxmodules,
            rng
    )
    #popsize is now 2
    @test populationtracker[1].Nvec[end] == length(populationtracker[1].cells)
    @test length(populationtracker[1].Nvec) === length(populationtracker[1].tvec) == 2
    #each cell has one mutation
    @test populationtracker[1].cells[1].data.mutations == populationtracker[1].cells[2].data.mutations == 1

    #check multilevel_simulation function
    rng = MersenneTwister(12)
    input = MultilevelBranchingInput(
            modulesize=4, 
            fixedmu=true, 
            b=0.1, 
            d=0.01,
            bdrate=0.01, 
            clonalmutations=0, 
            maxtime=20*365, 
            maxmodules=10,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = multilevel_simulation(SimpleTreeCell, input, rng)
    @test length(population) == 10
    population = multilevel_simulation(TreeCell, input, rng)
    @test length(population) == 10
    input = MultilevelBranchingInput(
            modulesize=4, 
            fixedmu=true, 
            b=1, 
            d=0,
            bdrate=0.1, 
            clonalmutations=0, 
            maxtime=365*100, 
            maxmodules=10,
            branchrate=3/365, 
            branchfraction=0.2, 
            μ=1
    )
    population = multilevel_simulation(SimpleTreeCell, input, rng)
    @test age(population) <= 365*50
    population = multilevel_simulation(TreeCell, input, rng)



end
