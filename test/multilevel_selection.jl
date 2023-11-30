@testset "cell" begin
    mt1 = SomaticEvolution.CellModule(
        Union{Cell, Nothing}[Cell([3, 4, 20], 2, 0), Cell([3, 5, 6, 8, 10, 17, 18, 19], 1, 0), Cell([3, 4, 21, 22, 23], 2, 0), Cell([1, 24], 1, 0)], 
        242.1,
        [0.0, 187.12],
        1, 0, WellMixed())
    mt2 = SomaticEvolution.CellModule(
        Union{Cell, Nothing}[Cell([1, 13, 25, 26, 27], 1, 0), Cell([1, 13], 1, 0), Cell([1, 13, 25, 28, 29], 1, 0)], 
        257.2, 
        [257.2],
        2, 1, WellMixed())
    mt3 = SomaticEvolution.CellModule(
        Union{Cell, Nothing}[Cell([1], 3, 0)],
        257.2,
        [257.2], 
        3, 2, WellMixed())

    @testset "updates" begin
        population = Population(
            [deepcopy(mt1)], 
            [deepcopy(mt2), deepcopy(mt3)], 
            [
                SomaticEvolution.Subclone(1, 0, 0.0, 5, 1.0, 0.1, 1.0, 0.1)
                SomaticEvolution.Subclone(2, 1, 240.0, 2, 1.5, 0.1, 1.5, 0.15)
                SomaticEvolution.Subclone(3, 1, 250.0, 1, 2.0, 0.1, 2.0, 0.2)
            ],
        )
        nsubclones = 4
        nsubclonescurrent = 3
        branchrate = 0.01
        t = 300
        nextID = 30
        μ, mutationdist = 1.0, :poisson
        #[(subclone id, number cells), (subclone id, number cells), ...]
        #Homeostatic modules: [(1, 2), (2, 2)]
        #Growing modules: [(1, 3)], [(3, 1)]
        @test SomaticEvolution.allclonetypes(population.growing_modules) == [1, 1, 1, 3]
        @test SomaticEvolution.get_moduleid_cellid(population.growing_modules, 1) == (1, 1)
        @test SomaticEvolution.get_moduleid_cellid(population.growing_modules, 3) == (1, 3)
        @test SomaticEvolution.get_moduleid_cellid(population.growing_modules, 4) == (2, 1)

        moduleid, cellid = SomaticEvolution.choose_module_cell(population.growing_modules, 1, rng)
        @test moduleid == 1
        @test cellid <=3
        moduleid, cellid = SomaticEvolution.choose_module_cell(population.growing_modules, 3, rng)
        @test moduleid == 2
        @test cellid == 1

        transitionrates = SomaticEvolution.get_selection_transitionrates(population, 0.01, nsubclones)
        @test round.(transitionrates, digits=2) == Float64[
            2.0, 3.0, 0.0, 0.0, #moran
            0.2, 0.3, 0.0, 0.0, #asymmetric
            3.0, 0.0, 2.0, 0.0, #birth
            0.3, 0.0, 0.1, 0.0, #death
            0.01 #newmodule
        ]

        #check choosing cell/module
        @testset "choose cell" begin
            homeostaticmoduleid, parentcellid = 
                SomaticEvolution.choose_module_cell(population.homeostatic_modules, 2, rng)
            homeostaticmodule = population.homeostatic_modules[homeostaticmoduleid]
            @test homeostaticmoduleid == 1
            @test homeostaticmodule.cells[parentcellid].clonetype == 2
            #choose cell to die from same module 
            deadcellid = 
                SomaticEvolution.choose_moran_deadcell(4, parentcellid, false, rng)
            @test  deadcellid != parentcellid
            growingmoduleid, parentcellid = 
                SomaticEvolution.choose_module_cell(population.growing_modules, 3, rng)
            @test growingmoduleid == 2
            @test SomaticEvolution.choose_homeostaticmodule(population, rng) == 1
        end

        @testset "transition rate updating" begin
        rng = MersenneTwister(12)
        #1. moran update
            parentcellid = 1
            deadcellid = 2
            homeostaticmodule = population.homeostatic_modules[1]
            parentcell = homeostaticmodule[parentcellid]
            deadcell = homeostaticmodule[deadcellid]    

            #implement cell division
            homesotaticmodule, subclones, nextID = 
                SomaticEvolution.celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, 
                    μ, mutationdist, rng)
            SomaticEvolution.celldeath!(homeostaticmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)

            transitionrates = 
                SomaticEvolution.update_selection_transitionrates_after_moran!(
                    transitionrates, population, getclonetype(parentcell), getclonetype(deadcell), nsubclones)   
            #Homeostatic modules: [(1, 1), (2, 3)]
            #Growing modules: [(1, 3)], [(3, 1)]
            @test length.(population.subclones) == [4, 3, 1]
            @test round.(transitionrates, digits=2) == Float64[
                1.0, 4.5, 0.0, 0.0, #moran
                0.1, 0.45, 0.0, 0.0, #asymmetric
                3.0, 0.0, 2.0, 0.0, #birth
                0.3, 0.0, 0.1, 0.0, #death
                0.01 #newmodule
            ]
            @test population.subclones[2].size == 3
            @test population.subclones[1].size == 4

        #2. birth update with transition to homeostasis
            parentcellid = 1
            transitiontohomeostasis = true
            growingmodule = population.growing_modules[1]
            parentcell = growingmodule[parentcellid]
            SomaticEvolution.move_module_to_homeostasis!(population, 1)
            growingmodule, subclones, nextID = 
                SomaticEvolution.celldivision!(growingmodule, population.subclones, parentcellid, t, nextID, μ, 
                    mutationdist, rng)
            transitionrates = 
                SomaticEvolution.update_selection_transitionrates_after_birth!(
                    transitionrates, population, getclonetype(parentcell), 
                    nsubclones, transitiontohomeostasis, branchrate, growingmodule)   

            #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
            #Growing modules: [(3, 1)]
            @test round.(transitionrates, digits=2) == Float64[
                5.0, 4.5, 0.0, 0.0, #moran
                0.5, 0.45, 0.0, 0.0, #asymmetric
                0.0, 0.0, 2.0, 0.0, #birth
                0.0, 0.0, 0.1, 0.0, #death
                0.02 #newmodule
            ]
            @test map(x -> x.size, population.subclones) == [5, 3, 1]

        #3. birth update with no transition to homeostasis and with mutation
            parentcellid = 1
            transitiontohomeostasis = false
            growingmodule = population.growing_modules[1]
            parentcell = growingmodule[parentcellid]
            growingmodule, subclones, nextID = SomaticEvolution.celldivision!(growingmodule, population.subclones, 
                parentcellid, t, nextID, μ, mutationdist, rng)
            transitionrates = 
                SomaticEvolution.update_selection_transitionrates_after_birth!(
                    transitionrates, population, growingmodule.cells[parentcellid].clonetype, 
                    nsubclones, transitiontohomeostasis, branchrate, growingmodule)   
            #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
            #Growing modules: [(3, 2)]
            @test round.(transitionrates, digits=2) == Float64[
                5.0, 4.5, 0.0, 0.0, #moran
                0.5, 0.45, 0.0, 0.0, #asymmetric
                0.0, 0.0, 4.0, 0.0, #birth
                0.0, 0.0, 0.2, 0.0, #death
                0.02 #newmodule
            ]
            @test population.subclones[3].size == 2

            selection = SelectionPredefined(Float64[0.5, 1.0, 0.5], Float64[1.0, 2.0, 3.0])
            moduletype = :growing
            population, growingmodule, nsubclonescurrent, transitionrates =
                SomaticEvolution.update_cellmutation!(population, growingmodule, nsubclonescurrent, nsubclones,
                    transitionrates, parentcell, selection, t, moduletype, rng)
            @test nsubclonescurrent == 4
            @test nsubclones == nsubclonescurrent
            #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
            #Growing modules: [(3, 1), (4, 1)]
            @test round.(transitionrates, digits=2) == Float64[
                5.0, 4.5, 0.0, 0.0, #moran
                0.5, 0.45, 0.0, 0.0, #asymmetric
                0.0, 0.0, 2.0, 1.5, #birth
                0.0, 0.0, 0.1, 0.1, #death
                0.02 #newmodule
            ]
            @test population.subclones[3].size == 1
            @test population.subclones[4].size == 1
                
        #4. death update
            deadcellid = 2
            growingmodule = population.growing_modules[1]
            deadcell = growingmodule[deadcellid]
            SomaticEvolution.celldeath!(growingmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
            transitionrates = 
                SomaticEvolution.update_selection_transitionrates_after_death!(transitionrates, 
                    population, getclonetype(deadcell), nsubclones)  
            #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
            #Growing modules: [(4, 1)]
            @test round.(transitionrates, digits=2) == Float64[
                5.0, 4.5, 0.0, 0.0, #moran
                0.5, 0.45, 0.0, 0.0, #asymmetric
                0.0, 0.0, 0.0, 1.5, #birth
                0.0, 0.0, 0.0, 0.1, #death
                0.02 #newmodule
            ]

        #5. module branching update (:split)
            branchinitsize = 2
            modulebranching = :split
            parentmoduleid = 2
            parentmodule = population.homeostatic_modules[parentmoduleid]
            parentmodule, newmodule, nextID = 
                SomaticEvolution.newmoduleformation!(parentmodule, population.subclones, 4, branchinitsize, t, rng; 
                    modulebranching, nextID, μ, mutationdist)
            push!(population.growing_modules, newmodule)
            SomaticEvolution.move_module_to_growing!(population, parentmoduleid)
            transitionrates = SomaticEvolution.update_selection_transitionrates_after_newmodule!(
                transitionrates, population, parentmodule, newmodule, nsubclones, branchrate, 4)
            #Homeostatic modules: [(1, 1), (2, 3)]
            #Growing modules: [(4, 1)], [(1, 2)], [(1, 2)]
            @test round.(transitionrates, digits=2) == Float64[
                1.0, 4.5, 0.0, 0.0, #moran
                0.1, 0.45, 0.0, 0.0, #asymmetric
                4.0, 0.0, 0.0, 1.5, #birth
                0.4, 0.0, 0.0, 0.1, #death
                0.01 #newmodule
            ]

        #6. module branching update (:withoutreplacement_nomutations)
            rng = MersenneTwister(12)
            branchinitsize = 1
            modulebranching = :withoutreplacement_nomutations
            parentmoduleid = 1
            parentmodule = population.homeostatic_modules[parentmoduleid]
            parentmodule, newmodule, nextID = 
                SomaticEvolution.newmoduleformation!(parentmodule, population.subclones, 5, branchinitsize, t, rng; 
                    modulebranching, nextID, μ, mutationdist)
            push!(population.growing_modules, newmodule)
            # @show newmodule[1].clonetype == 1 (this is random)
            transitionrates = SomaticEvolution.update_selection_transitionrates_after_newmodule!(
                transitionrates, population, parentmodule, newmodule, nsubclones, branchrate, 4)
            #Homeostatic modules: [(1, 1), (2, 3)]
            #Growing modules: [(4, 1)], [(1, 2)], [(1, 2)], [(1, 1)]
            @test round.(transitionrates, digits=2) == Float64[
                1.0, 4.5, 0.0, 0.0, #moran
                0.1, 0.45, 0.0, 0.0, #asymmetric
                5.0, 0.0, 0.0, 1.5, #birth
                0.5, 0.0, 0.0, 0.1, #death
                0.01 #newmodule
            ]

        #7. kill growing module
        deadmoduleid = 3
        deadmodule = population[deadmoduleid]
        SomaticEvolution.moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
        @test length.(population.subclones) == [4, 3, 0, 1]
        transitionrates = SomaticEvolution.update_selection_transitionrates_after_moduledeath!(
            transitionrates, population, deadmodule, nsubclones, branchrate, 4)
        #Homeostatic modules: [(1, 1), (2, 3)]
        #Growing modules: [(4, 1)], [(1, 2)], [(1, 1)]
        @test round.(transitionrates, digits=2) == Float64[
            1.0, 4.5, 0.0, 0.0, #moran
            0.1, 0.45, 0.0, 0.0, #asymmetric
            3.0, 0.0, 0.0, 1.5, #birth
            0.3, 0.0, 0.0, 0.1, #death
            0.01 #newmodule
        ]
        end        
    end


    @testset "simulate full predefined selection" begin

        #check runsimulation function
        rng = MersenneTwister(12)
        input = MultilevelBranchingMoranInput(
                modulesize=4, 
                mutationdist=:fixed, 
                birthrate=10, 
                deathrate=0.01,
                moranrate=10, 
                clonalmutations=0, 
                tmax=3, 
                maxmodules=3,
                branchrate=3, 
                branchinitsize=1, 
                μ=1,
        )
        selection = SelectionPredefined([0.2, 1.0], [0.5, 2.0])
        population = runsimulation(Cell, input, selection, rng)
        @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
        @test length(population) == 3
        @test age(population) < 365 * 3
    end

    @testset "simulate full selection distribution" begin

        #check runsimulation function
        rng = MersenneTwister(12)
        input = MultilevelBranchingMoranInput(
                modulesize=4, 
                mutationdist=:fixed, 
                birthrate=10, 
                deathrate=0.01,
                moranrate=10, 
                clonalmutations=0, 
                tmax=3, 
                maxmodules=3,
                branchrate=3, 
                branchinitsize=1, 
                μ=1,
        )
        mean, std = 0.1, 0.02
        shape, scale = mean^2/std^2, std^2/mean
        selection = SelectionDistribution(
            Gamma(shape, scale),
            0.1,
            50
        )
        simulation = runsimulation(Cell, input, selection, rng)
        population = simulation.output
        @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
        @test length(population) == 3
        @test age(population) < 365 * 3
        @test (length.(population.subclones) == 
            counts(SomaticEvolution.allclonetypes(population), 1:length(population.subclones)))
    end



end

# @testset "simple tree" begin
#     rng = MersenneTwister(12)
#     population = SomaticEvolution.initialize_population(
#         SimpleTreeCell, WellMixed, 0, 3, 3, 100, 1, 100, 50, 2, rng=rng)

#     @testset "updates" begin
#         nsubclones = 3
#         nsubclonescurrent = 1
#         branchrate = 1
#         t = 1
#         nextID = 2
#         μ, mutationdist = 1.0, :poisson
#         #[(subclone id, number cells), (subclone id, number cells), ...]
#         #Homeostatic modules: [(1, 3), (1, 3)]
#         @test SomaticEvolution.allclonetypes(population.growing_modules) == [1, 1, 1, 1, 1, 1]
#         @test SomaticEvolution.get_moduleid_cellid(population.growing_modules, 1) == (1, 1)
#         @test SomaticEvolution.get_moduleid_cellid(population.growing_modules, 3) == (1, 3)
#         @test SomaticEvolution.get_moduleid_cellid(population.growing_modules, 4) == (2, 1)

#         moduleid, cellid = SomaticEvolution.choose_module_cell(population.growing_modules, 1, rng)
#         @test moduleid == 1
#         @test cellid <=3
#         moduleid, cellid = SomaticEvolution.choose_module_cell(population.growing_modules, 3, rng)
#         @test moduleid == 2
#         @test cellid == 1

#         transitionrates = SomaticEvolution.get_selection_transitionrates(population, 0.01, nsubclones)
#         @test round.(transitionrates, digits=2) == Float64[
#             600.0, 0.0, 0.0, #moran
#             300.0, 0.3, 0.0, #asymmetric
#             600.0, 0.0, 2.0, #birth
#             6.0, 0.0, 0.1, #death
#             2.0 #newmodule
#         ]

#         @testset "transition rate updating" begin
#         rng = MersenneTwister(12)
#         #1. moran update
#             parentcellid = 1
#             deadcellid = 2
#             homeostaticmodule = population.homeostatic_modules[1]
#             parentcell = homeostaticmodule[parentcellid]
#             deadcell = homeostaticmodule[deadcellid]    

#             #implement cell division
#             homesotaticmodule, subclones, nextID = 
#                 SomaticEvolution.celldivision!(homeostaticmodule, population.subclones, parentcellid, t, nextID, 
#                     μ, mutationdist, rng)
#             SomaticEvolution.celldeath!(homeostaticmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)

#             transitionrates = 
#                 SomaticEvolution.update_selection_transitionrates_after_moran!(
#                     transitionrates, population, getclonetype(parentcell), getclonetype(deadcell), nsubclones)   
#             #Homeostatic modules: [(1, 1), (2, 3)]
#             #Growing modules:
#             @test round.(transitionrates, digits=2) == Float64[
#                 1.0, 4.5, 0.0, 0.0, #moran
#                 0.1, 0.45, 0.0, 0.0, #asymmetric
#                 3.0, 0.0, 2.0, 0.0, #birth
#                 0.3, 0.0, 0.1, 0.0, #death
#                 0.01 #newmodule
#             ]
#             @test population.subclones[2].size == 3
#             @test population.subclones[1].size == 4

#         #2. birth update with transition to homeostasis
#             parentcellid = 1
#             transitiontohomeostasis = true
#             growingmodule = population.growing_modules[1]
#             parentcell = growingmodule[parentcellid]
#             SomaticEvolution.move_module_to_homeostasis!(population, 1)
#             growingmodule, subclones, nextID = 
#                 SomaticEvolution.celldivision!(growingmodule, population.subclones, parentcellid, t, nextID, μ, 
#                     mutationdist, rng)
#             transitionrates = 
#                 SomaticEvolution.update_selection_transitionrates_after_birth!(
#                     transitionrates, population, getclonetype(parentcell), 
#                     nsubclones, transitiontohomeostasis, branchrate, growingmodule)   

#             #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
#             #Growing modules: [(3, 1)]
#             @test round.(transitionrates, digits=2) == Float64[
#                 5.0, 4.5, 0.0, 0.0, #moran
#                 0.5, 0.45, 0.0, 0.0, #asymmetric
#                 0.0, 0.0, 2.0, 0.0, #birth
#                 0.0, 0.0, 0.1, 0.0, #death
#                 0.02 #newmodule
#             ]
#             @test map(x -> x.size, population.subclones) == [5, 3, 1]

#         #3. birth update with no transition to homeostasis and with mutation
#             parentcellid = 1
#             transitiontohomeostasis = false
#             growingmodule = population.growing_modules[1]
#             parentcell = growingmodule[parentcellid]
#             growingmodule, subclones, nextID = SomaticEvolution.celldivision!(growingmodule, population.subclones, 
#                 parentcellid, t, nextID, μ, mutationdist, rng)
#             transitionrates = 
#                 SomaticEvolution.update_selection_transitionrates_after_birth!(
#                     transitionrates, population, growingmodule.cells[parentcellid].clonetype, 
#                     nsubclones, transitiontohomeostasis, branchrate, growingmodule)   
#             #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
#             #Growing modules: [(3, 2)]
#             @test round.(transitionrates, digits=2) == Float64[
#                 5.0, 4.5, 0.0, 0.0, #moran
#                 0.5, 0.45, 0.0, 0.0, #asymmetric
#                 0.0, 0.0, 4.0, 0.0, #birth
#                 0.0, 0.0, 0.2, 0.0, #death
#                 0.02 #newmodule
#             ]
#             @test population.subclones[3].size == 2

#             moduletype = :growing
#             population, growingmodule, nsubclonescurrent, transitionrates =
#                 SomaticEvolution.update_cellmutation!(population, growingmodule, nsubclonescurrent, nsubclones,
#                     transitionrates, parentcell, 0.5, t, moduletype)
#             @test nsubclonescurrent == 4
#             @test nsubclones == nsubclonescurrent
#             #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
#             #Growing modules: [(3, 1), (4, 1)]
#             @test round.(transitionrates, digits=2) == Float64[
#                 5.0, 4.5, 0.0, 0.0, #moran
#                 0.5, 0.45, 0.0, 0.0, #asymmetric
#                 0.0, 0.0, 2.0, 1.5, #birth
#                 0.0, 0.0, 0.1, 0.1, #death
#                 0.02 #newmodule
#             ]
#             @test population.subclones[3].size == 1
#             @test population.subclones[4].size == 1
                
#         #4. death update
#             deadcellid = 2
#             growingmodule = population.growing_modules[1]
#             deadcell = growingmodule[deadcellid]
#             SomaticEvolution.celldeath!(growingmodule, population.subclones, deadcellid, t, μ, mutationdist, rng)
#             transitionrates = 
#                 SomaticEvolution.update_selection_transitionrates_after_death!(transitionrates, 
#                     population, getclonetype(deadcell), nsubclones)  
#             #Homeostatic modules: [(1, 1), (2, 3)], [(1, 4)]
#             #Growing modules: [(4, 1)]
#             @test round.(transitionrates, digits=2) == Float64[
#                 5.0, 4.5, 0.0, 0.0, #moran
#                 0.5, 0.45, 0.0, 0.0, #asymmetric
#                 0.0, 0.0, 0.0, 1.5, #birth
#                 0.0, 0.0, 0.0, 0.1, #death
#                 0.02 #newmodule
#             ]

#         #5. module branching update (:split)
#             branchinitsize = 2
#             modulebranching = :split
#             parentmoduleid = 2
#             parentmodule = population.homeostatic_modules[parentmoduleid]
#             parentmodule, newmodule, nextID = 
#                 SomaticEvolution.newmoduleformation!(parentmodule, population.subclones, 4, branchinitsize, t, rng; 
#                     modulebranching, nextID, μ, mutationdist)
#             push!(population.growing_modules, newmodule)
#             SomaticEvolution.move_module_to_growing!(population, parentmoduleid)
#             transitionrates = SomaticEvolution.update_selection_transitionrates_after_newmodule!(
#                 transitionrates, population, parentmodule, newmodule, nsubclones, branchrate, 4)
#             #Homeostatic modules: [(1, 1), (2, 3)]
#             #Growing modules: [(4, 1)], [(1, 2)], [(1, 2)]
#             @test round.(transitionrates, digits=2) == Float64[
#                 1.0, 4.5, 0.0, 0.0, #moran
#                 0.1, 0.45, 0.0, 0.0, #asymmetric
#                 4.0, 0.0, 0.0, 1.5, #birth
#                 0.4, 0.0, 0.0, 0.1, #death
#                 0.01 #newmodule
#             ]

#         #6. module branching update (:withoutreplacement_nomutations)
#             rng = MersenneTwister(12)
#             branchinitsize = 1
#             modulebranching = :withoutreplacement_nomutations
#             parentmoduleid = 1
#             parentmodule = population.homeostatic_modules[parentmoduleid]
#             parentmodule, newmodule, nextID = 
#                 SomaticEvolution.newmoduleformation!(parentmodule, population.subclones, 5, branchinitsize, t, rng; 
#                     modulebranching, nextID, μ, mutationdist)
#             push!(population.growing_modules, newmodule)
#             # @show newmodule[1].clonetype == 1 (this is random)
#             transitionrates = SomaticEvolution.update_selection_transitionrates_after_newmodule!(
#                 transitionrates, population, parentmodule, newmodule, nsubclones, branchrate, 4)
#             #Homeostatic modules: [(1, 1), (2, 3)]
#             #Growing modules: [(4, 1)], [(1, 2)], [(1, 2)], [(1, 1)]
#             @test round.(transitionrates, digits=2) == Float64[
#                 1.0, 4.5, 0.0, 0.0, #moran
#                 0.1, 0.45, 0.0, 0.0, #asymmetric
#                 5.0, 0.0, 0.0, 1.5, #birth
#                 0.5, 0.0, 0.0, 0.1, #death
#                 0.01 #newmodule
#             ]

#         #7. kill growing module
#         deadmoduleid = 3
#         deadmodule = population[deadmoduleid]
#         SomaticEvolution.moduledeath!(population, deadmoduleid, t, μ, mutationdist, rng)
#         transitionrates = SomaticEvolution.update_selection_transitionrates_after_moduledeath!(
#             transitionrates, population, deadmodule, nsubclones, branchrate, 4)
#         #Homeostatic modules: [(1, 1), (2, 3)]
#         #Growing modules: [(4, 1)], [(1, 2)], [(1, 1)]
#         @test round.(transitionrates, digits=2) == Float64[
#             1.0, 4.5, 0.0, 0.0, #moran
#             0.1, 0.45, 0.0, 0.0, #asymmetric
#             3.0, 0.0, 0.0, 1.5, #birth
#             0.3, 0.0, 0.0, 0.1, #death
#             0.01 #newmodule
#         ]
#         end        
#     end


#     # @testset "simulate full" begin

#     #     #check runsimulation function
#     #     rng = MersenneTwister(12)
#     #     input = MultilevelBranchingMoranInput(
#     #             modulesize=4, 
#     #             mutationdist=:fixed, 
#     #             birthrate=10, 
#     #             deathrate=0.01,
#     #             moranrate=10, 
#     #             clonalmutations=0, 
#     #             tmax=3, 
#     #             maxmodules=3,
#     #             branchrate=3, 
#     #             branchinitsize=1, 
#     #             μ=1,
#     #             mutant_selection=[0.2, 1.0],
#     #             mutant_time=[0.5, 2.0]
#     #     )
#     #     population = runsimulation(Cell, input, rng)
#     #     @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
#     #     @test length(population) == 3
#     #     @test age(population) < 365 * 3

#     #     input = MultilevelBranchingMoranInput(
#     #         modulesize=4, 
#     #         mutationdist=:fixed, 
#     #         birthrate=10, 
#     #         deathrate=0.01,
#     #         moranrate=10, 
#     #         clonalmutations=0, 
#     #         tmax=3, 
#     #         maxmodules=3,
#     #         branchrate=3, 
#     #         branchinitsize=1, 
#     #         μ=1,
#     #         mutant_selection=[5.0],
#     #         mutant_time=[0.5]
#     # )
#     # population = runsimulation(SimpleTreeCell, input, rng)
#     # @test sort(moduleid.(population)) == sort(unique(moduleid.(population)))
#     # @test length(population) == 3
#     # @test age(population) < 365 * 3
#     # end

# end