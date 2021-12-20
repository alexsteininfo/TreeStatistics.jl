rng  = MersenneTwister(12)

inputs = InputParameters{MoranInput}[]
for N in [100, 1000]     
    for numclones in [0, 1, 2]
        for selection in ([0.1,0.4],[1.0,3.0],[2.0,10.0])
            for tevent in ([1, 10], [10, 15])
                IP = InputParameters{MoranInput}(N=100, numclones=numclones, μ=100, fixedmu=false, 
                    clonalmutations=0, tmax=20, bdrate=1, selection=selection[1:numclones],
                    tevent = tevent[1:numclones])
                push!(inputs, IP)
            end
        end
    end
end

for IP in inputs
    simtracker = SomaticEvolution.initializesim(IP.siminput)
    @test sum(simtracker.clonesize) == IP.siminput.N
    simtracker = 
        SomaticEvolution.moranprocess!(simtracker, IP.siminput.bdrate, IP.siminput.tmax, 
            IP.siminput.μ, rng, numclones = IP.siminput.numclones, fixedmu = IP.siminput.fixedmu, 
            selection = IP.siminput.selection, tevent = IP.siminput.tevent)
    @test length(simtracker.subclones) == IP.siminput.numclones
    simtracker, simresults = 
        SomaticEvolution.processresults!(simtracker, IP.siminput.N, IP.siminput.numclones, IP.siminput.μ, 
                    IP.siminput.fixedmu, IP.siminput.clonalmutations, 
                    IP.ploidy, 0, 1, rng)
    @test length(simresults.subclones) == IP.siminput.numclones
    for i in 1:IP.siminput.numclones
        @test (sum([cell.clonetype for cell in simresults.cells] .== i+1)/IP.siminput.N
                    == simresults.subclones[i].freq)
    end                
                    

    
end