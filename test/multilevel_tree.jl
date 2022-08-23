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

module1 = SomaticEvolution.TreeModule(Int64[], Float64[], [node19, node24, node25], SomaticEvolution.CloneTracker[], 1, 0)
module2 = SomaticEvolution.TreeModule(Int64[], Float64[], [node26, node27, node28], SomaticEvolution.CloneTracker[], 2, 1)
module3 = SomaticEvolution.TreeModule(Int64[], Float64[], [node29, node30, node23], SomaticEvolution.CloneTracker[], 3, 2)

population = [module1, module2, module3]
@testset "pairwise differences" begin
    @test SomaticEvolution.findMRCA(module1) == mrca1
    @test SomaticEvolution.findMRCA(module2) == mrca2
    @test SomaticEvolution.findMRCA(module3) == mrca3
    @test pairwise_fixed_differences_clonal(population) == (Dict(68 => 1, 102 => 1, 54 => 1), Dict(42 => 1, 32 => 1, 66 => 1))
end












