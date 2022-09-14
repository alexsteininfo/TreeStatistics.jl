@testset "neutral runsimulation" begin
    rng = MersenneTwister(100)
    input = BranchingInput(
        Nmax=10, 
        mutationdist=:poisson, 
        b=1, 
        d=0.0,
        clonalmutations=0, 
        numclones=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test length(simulation.output) == 10
    @test simulation.output.Nvec[end] == length(simulation.output.cells)

    tmax=10
    input = MoranInput(
        N=10, 
        tmax=tmax,
        mutationdist=:poisson, 
        bdrate=1.0,
        clonalmutations=0, 
        numclones=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test all(simulation.output.Nvec .== 10)
    @test simulation.output.tvec[end] <=tmax
    tmax=10
    input = BranchingMoranInput(
        Nmax=10, 
        tmax=tmax,
        mutationdist=:poisson, 
        b=1, 
        d=0.0,
        clonalmutations=0, 
        numclones=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test simulation.output.Nvec[end]== 10
    @test simulation.output.tvec[end] <= tmax
end

@testset "cell transitions" begin
    
end

subclones = [
    SomaticEvolution.CloneTracker(
        1, 1, 1.4437821887595856, [2], 2, 1, 1.0, 9
    ),
    SomaticEvolution.CloneTracker(
        1, 1, 1.7344043713478965, [1, 4], 3, 2, 1.6666666666666667, 90
    )
]
N = 100
@test SomaticEvolution.getclonesize(N, subclones) == [1, 9, 90]

#check birth and death rates are calculated correctly
b, d, selection = 1, 0, [1, 2]
brates, drates = SomaticEvolution.set_branching_birthdeath_rates(b, d, selection) 
@test brates == [1, 2, 3]
@test drates == [0, 0, 0]