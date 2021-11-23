"""
    show(sresult::Simulation)

Print out summary of simulation.
"""
function Base.show(io::IO, sresult::Simulation{BranchingInput})
    @printf("===================================================================\n")
    @printf("Branching process \n")
    @printf("Input parameters: \n")
    @printf("\t Mutation rate: %.2f\n", sresult.input.siminput.μ)
    @printf("\t Death rate of host population: %.2f\n", sresult.input.siminput.d)
    @printf("\t Effective mutation rate (μ/β): %.2f\n",(sresult.input.siminput.μ / 
        ((sresult.input.siminput.b-sresult.input.siminput.d)/sresult.input.siminput.b)))
    @printf("\t Number of clonal mutation: %d\n", sresult.input.siminput.clonalmutations)
    @printf("\t Number of subclones: %d\n\n", sresult.input.siminput.numclones)
    if sresult.input.siminput.numclones > 0
        for i in 1:length(sresult.output.clonefreq)
            @printf("Subclone %d \n", i)
            @printf("\tFrequency: %.2f\n", sresult.output.clonefreq[i])
            @printf("\tNumber of mutations in subclone: %d\n", sresult.output.subclonalmutations[i])
            @printf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i])
            @printf("\tTime subclone emerges: %.2f\n", sresult.output.clonetime[i])
            @printf("\tNumber of divisions: %d\n", sresult.output.Ndivisions[i])
            @printf("\tAverage number of divisions per cell: %.2f\n", sresult.output.avdivisions[i])
            @printf("\tPopulation size when subclone emerges: %d\n", sresult.output.cloneN[i])
            @printf("\tParent of subclone (0 is host): %d\n\n", sresult.output.clonetype[i] - 1)
        end
    else
        @printf("No clones, tumour growth was neutral\n\n")
    end

end

"""
    show(sresult::Simulation)

Print out summary of simulation.
"""
function Base.show(io::IO, sresult::Simulation{MoranInput})
    @printf("===================================================================\n")
    @printf("Moran process \n")
    @printf("Input parameters: \n")
    @printf("\t Mutation rate: %.2f\n", sresult.input.siminput.μ)
    @printf("\t Birth/death rate of host population: %.2f\n", sresult.input.siminput.bdrate)
    @printf("\t Number of clonal mutation: %d\n", sresult.input.siminput.clonalmutations)
    @printf("\t Number of subclones: %d\n\n", sresult.input.siminput.numclones)
    if sresult.input.siminput.numclones > 0
        for i in 1:length(sresult.output.clonefreq)
            @printf("Subclone %d \n", i)
            @printf("\tFrequency: %.2f\n", sresult.output.clonefreq[i])
            @printf("\tNumber of mutations in subclone: %d\n", sresult.output.subclonalmutations[i])
            @printf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i])
            @printf("\tTime subclone emerges: %.2f\n", sresult.output.clonetime[i])
            @printf("\tNumber of divisions: %d\n", sresult.output.Ndivisions[i])
            @printf("\tAverage number of divisions per cell: %.2f\n", sresult.output.avdivisions[i])
            @printf("\tPopulation size when subclone emerges: %d\n", sresult.output.cloneN[i])
            @printf("\tParent of subclone (0 is host): %d\n\n", sresult.output.clonetype[i] - 1)
        end
    else
        @printf("No clones, tumour growth was neutral\n\n")
    end
end

function Base.show(io::IO, sresult::MultiSimulation)
    @printf("===================================================================\n")
    @printf("Number of simulations = %d\n\n", length(sresult.output))
    @printf("Input parameters: \n")
    @printf("\t Mutation rate: %.2f\n", sresult.input.siminput.μ)
    @printf("\t Death rate of host population: %.2f\n", sresult.input.siminput.d)
    @printf("\t Effective mutation rate (μ/β): %.2f\n",(sresult.input.siminput.μ / 
      ((sresult.input.siminput.b-sresult.input.siminput.d)/sresult.input.siminput.b)))
    @printf("\t Number of clonal mutation: %d\n", sresult.input.siminput.clonalmutations)
    @printf("\t Number of subclones: %d\n\n", sresult.input.siminput.numclones)  

    if sresult.input.siminput.numclones > 0
        for i in 1:sresult.input.siminput.numclones
            @printf("Subclone %d \n", i)
            @printf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i])
            @printf("\tTime subclone emerges (population doublings): %.2f\n", 
                sresult.input.siminput.tevent[i])
            @printf("\tAverage frequency: %.2f (std: %.2f)\n", 
                mean([simoutput.clonefreq[i] for simoutput in sresult.output]),
                std([simoutput.clonefreq[i] for simoutput in sresult.output]))
            @printf("\tAverage number of mutations in subclone: %d (std: %.2f)\n", 
                mean([simoutput.subclonalmutations[i] for simoutput in sresult.output]),
                std([simoutput.subclonalmutations[i] for simoutput in sresult.output]))
            @printf("\tAverage population size when subclone emerges: %d (std: %.2f)\n", 
                mean([simoutput.cloneN[i] for simoutput in sresult.output]),
                std([simoutput.cloneN[i] for simoutput in sresult.output]))
            parenttypes = [simoutput.clonetype[i] - 1 for simoutput in sresult.output]
            parenttypes = sort(collect(countmap(parenttypes)))
            @printf("\tDistribution of subclone parents (0 is host): ")
            s = join([@sprintf("%d => %d", parent, freq) for (parent, freq) in parenttypes]
                        , ", ")
            @printf("%s\n\n",s)
        end
    else
        @printf("No clones, tumour growth was neutral\n\n")
    end

end

function Base.show(io::IO, sresultlist::Vector{MultiSimulation{T}} where T <: SimulationInput)
    @printf("%d sets of simulations\n", length(sresultlist))
    for sresult in sresultlist
        Base.show(io, sresult)
    end
end

function saveinput(sresult::Simulation{BranchingInput}, filename)
    open(filename, "w") do io
        write(io, @sprintf("Simulation => branching process\n\n"))
        write(io, @sprintf("Input parameters:\n"))
        write(io, @sprintf("\tDetection limit = %.3f\n", sresult.input.detectionlimit))
        write(io, @sprintf("\tPloidy = %d\n",sresult.input.ploidy))
        write(io, @sprintf("\tRead depth = %d\n",sresult.input.read_depth))
        write(io, @sprintf("\tρ = %.2f\n",sresult.input.ρ))
        write(io, @sprintf("\tCellularity = %.2f\n\n",sresult.input.cellularity))
        write(io, @sprintf("\tMax population size = %d\n",sresult.input.siminput.Nmax))
        write(io, @sprintf("\tMutation rate (per division) = %d\n",sresult.input.siminput.μ))
        write(io, @sprintf("\tClonal mutations = %d\n",sresult.input.siminput.clonalmutations))
        write(io, @sprintf("\tBirth rate (wild type)= %.2f\n",sresult.input.siminput.b))
        write(io, @sprintf("\tDeath rate (wild type)= %.2f\n",sresult.input.siminput.d))
        write(io, @sprintf("\tMax size before mutations stop accumulating (will not be detected) = %.2f\n\n",sresult.input.siminput.maxclonesize))
        write(io, @sprintf("\tNumber of subclones = %d\n",sresult.input.siminput.numclones))
        if sresult.input.siminput.numclones > 0
            ss = join([@sprintf("%.2f", s) for s in sresult.input.siminput.selection], ", ")
            write(io, @sprintf("\tSelection strength = %s\n", ss))
            ts = join([@sprintf("%.2f", t) for t in sresult.input.siminput.tevent], ", ")
            write(io, @sprintf("\tSubclone event time = %s\n\n", ts))
            write(io, @sprintf("Subclone output:\n"))
            for i in 1:length(sresult.output.clonefreq)
                write(io, @sprintf("Subclone %d: \n", i))
                write(io, @sprintf("\tFrequency: %.2f\n", sresult.output.clonefreq[i]))
                write(io, @sprintf("\tNumber of mutations in subclone: %d\n", sresult.output.subclonalmutations[i]))
                write(io, @sprintf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i]))
                write(io, @sprintf("\tTime subclone emerges (population doublings): %.2f\n", sresult.output.clonetime[i]))
                write(io, @sprintf("\tNumber of divisions: %d\n", sresult.output.Ndivisions[i]))
                write(io, @sprintf("\tAverage number of divisions per cell: %.2f\n", sresult.output.avdivisions[i]))
                write(io, @sprintf("\tPopulation size when subclone emerges: %d\n", sresult.output.cloneN[i]))
                write(io, @sprintf("\tParent of subclone (0 is host): %d\n\n", sresult.output.clonetype[i] - 1))
            end
        else
            @printf("No clones, tumour growth was neutral\n\n")
        end
        

    end

end
    