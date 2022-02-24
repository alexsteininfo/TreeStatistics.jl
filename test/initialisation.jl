@testset "population initialization" begin
    for clonalmutations in [0,100]
        pop1 = SomaticEvolution.initialize_population(100, clonalmutations=clonalmutations)
        pop2 = SomaticEvolution.initialize_population(clonalmutations=clonalmutations)
        for pop in (pop1, pop2)
            @test length(pop) == 1
            @test pop[1].Nvec[1] == 1
            @test length(pop[1].Nvec) == length(pop[1].tvec) == 1
            @test pop[1].tvec[1] == 0.0
            @test length(pop[1].cells) == 1
            @test length(pop[1].cells[1].mutations) == clonalmutations
            @test pop[1].cells[1].clonetype == 1
            @test length(pop[1].subclones) == 0
            @test pop[1].id == 1
        end
    end
end

@testset "module initialization" begin
    @testset "branching" begin
        for clonalmutations in [0,100]
            module1 = SomaticEvolution.initializesim_branching(100, clonalmutations=clonalmutations)
            module2 = SomaticEvolution.initializesim_branching(clonalmutations=clonalmutations)
            for initmodule in (module1, module2)
                @test initmodule.Nvec[1] == 1
                @test length(initmodule.Nvec) == length(initmodule.tvec) == 1
                @test initmodule.tvec[1] == 0.0
                @test length(initmodule.cells) == 1
                @test length(initmodule.cells[1].mutations) == clonalmutations
                @test initmodule.cells[1].clonetype == 1
                @test length(initmodule.subclones) == 0
                @test initmodule.id == 1
            end
        end

        @testset "moran" begin
            for clonalmutations in [0,100]
                N = 100
                initmodule = SomaticEvolution.initializesim_moran(N, clonalmutations=clonalmutations)
                @test initmodule.Nvec[1] == N
                @test length(initmodule.Nvec) == length(initmodule.tvec) == 1
                @test initmodule.tvec[1] == 0.0
                @test length(initmodule.cells) == N
                allmuts = [cell.mutations for cell in initmodule.cells]
                @test all(x -> x == allmuts[1], allmuts)
                @test length(initmodule.cells[1].mutations) == clonalmutations
                @test initmodule.cells[1].clonetype == 1
                @test length(initmodule.subclones) == 0
                @test initmodule.id == 1
            end
        end

        @testset "sampled" begin
            subclones = SomaticEvolution.CloneTracker[]
            cells = Cell[]
            push!(cells, Cell([1,2,3,4], 1))
            push!(cells, Cell([1,2,3], 1))
            id = 1
            parentid = 0
            inittime = 5
            initmodule = SomaticEvolution.initializesim_from_cells(cells, subclones, id,
                parentid, inittime=inittime)
            @test initmodule.Nvec[1] == length(cells)
            @test length(initmodule.Nvec) == length(initmodule.tvec) == 1
            @test initmodule.tvec[1] == inittime
            @test length(initmodule.cells) == length(cells)
            @test initmodule.cells[1].mutations == [1,2,3,4]
            @test initmodule.cells[1].clonetype == initmodule.cells[2].clonetype == 1
            @test length(initmodule.subclones) == length(subclones)
            @test initmodule.id == id
        end
    end
end
