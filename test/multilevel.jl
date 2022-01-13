Nmax = 100
IP = InputParameters{BranchingMoranInput}(Nmax=Nmax, numclones=0, μ=100, fixedmu=false, bdrate=log(2),
    clonalmutations=0, tmax=10)

moduletracker = SomaticEvolution.initialize_population(
    IP.siminput.Nmax, 
    clonalmutations=IP.siminput.clonalmutations
    )[1]

branchinitsize = 20
modulebranchingrate = 1/5
maxage = 50
rng = MersenneTwister(12)

moduletracker, newmoduletracker =
    SomaticEvolution.module_simulate_to_branching(
        moduletracker, 
        IP,  
        branchinitsize,
        2,
        modulebranchingrate,
        maxage,
        rng
    )

@testset "simulate till branching" begin
    @test moduletracker.tvec[end] <= maxage
    @test moduletracker.Nvec[end] == Nmax - branchinitsize
    @test newmoduletracker.Nvec[end] == branchinitsize
    @test newmoduletracker.tvec[end] == moduletracker.tvec[end]
    
end