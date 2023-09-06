@testset "module initialization" begin
    @testset "branching" begin
        for clonalmutations in [0,100]
            module1 = SomaticEvolution.initialize(Cell, WellMixed, clonalmutations, 1)
            module2 = SomaticEvolution.initialize(Cell, WellMixed, clonalmutations, 1)
            for initmodule in (module1, module2)
                @test length(initmodule) == 1
                @test initmodule.t == 0.0
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
                initmodule = SomaticEvolution.initialize(Cell, WellMixed, clonalmutations, N)
                @test length(initmodule) == N
                @test initmodule.branchtimes[1] == 0.0
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
            push!(cells, Cell([1,2,3,4], 1, 0.2, 1, 0))
            push!(cells, Cell([1,2,3], 1, 0.2, 2, 1))
            id = 1
            parentid = 0
            inittime = 5
            initmodule = SomaticEvolution.new_module_from_cells(cells, inittime, [inittime], subclones, id,
                parentid, WellMixed())
            @test initmodule.branchtimes[1] == inittime
            @test length(initmodule.cells) == length(cells)
            @test initmodule.cells[1].mutations == [1,2,3,4]
            @test initmodule.cells[1].clonetype == initmodule.cells[2].clonetype == 1
            @test length(initmodule.subclones) == length(subclones)
            @test initmodule.id == id
        end
    end
end
