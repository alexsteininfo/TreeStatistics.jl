"""
    show(sresult::Simulation)

Print out summary of simulation.
"""
function Base.show(io::IO, sresult::Simulation)
    @printf("===================================================================\n")
    @printf("Branching process \n")
    @printf("Input parameters: \n")
    @printf("\t Mutation rate: %.2f\n", sresult.input.siminput.μ)
    _showrates(io, sresult)
    @printf("\t Number of clonal mutations: %d\n", sresult.input.siminput.clonalmutations)
    @printf("\t Number of subclones: %d\n\n", sresult.input.siminput.numclones)
    if sresult.input.siminput.numclones > 0
        for (i, subclone) in enumerate(sresult.output.subclones)
            @printf("Subclone %d \n", i)
            @printf("\tFrequency: %.2f\n", subclone.freq)
            @printf("\tNumber of mutations in subclone: %d\n", subclone.mutations)
            @printf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i])
            @printf("\tTime subclone emerges: %.2f\n", subclone.time)
            @printf("\tNumber of divisions: %d\n", subclone.Ndivisions)
            @printf("\tAverage number of divisions per cell: %.2f\n", subclone.avdivisions)
            @printf("\tPopulation size when subclone emerges: %d\n", subclone.N0)
            @printf("\tParent of subclone (0 is host): %d\n\n", subclone.parenttype-1)
        end
    else
        @printf("No clones, tumour growth was neutral\n\n")
    end

end

function _showrates(io::IO, sresult::Simulation{BranchingInput})
    @printf("\t Birth rate of host population: %.2f\n", sresult.input.siminput.b)
    @printf("\t Death rate of host population: %.2f\n", sresult.input.siminput.d)
    @printf("\t Effective mutation rate (μ/β): %.2f\n",(sresult.input.siminput.μ / 
        ((sresult.input.siminput.b-sresult.input.siminput.d)/sresult.input.siminput.b)))
end

function _showrates(io::IO, sresult::Simulation{MoranInput})
    @printf("\t Birth/death rate of host population: %.2f\n", sresult.input.siminput.bdrate)
end
"""
    show(sresult::Simulation)

Print out summary of simulation.
"""
# function Base.show(io::IO, sresult::Simulation{MoranInput})
#     @printf("===================================================================\n")
#     @printf("Moran process \n")
#     @printf("Input parameters: \n")
#     @printf("\t Mutation rate: %.2f\n", sresult.input.siminput.μ)
#     @printf("\t Birth/death rate of host population: %.2f\n", sresult.input.siminput.bdrate)
#     @printf("\t Number of clonal mutation: %d\n", sresult.input.siminput.clonalmutations)
#     @printf("\t Number of subclones: %d\n\n", sresult.input.siminput.numclones)
#     if sresult.input.siminput.numclones > 0
#         for i in 1:length(sresult.output.clonefreq)
#             @printf("Subclone %d \n", i)
#             @printf("\tFrequency: %.2f\n", sresult.output.clonefreq[i])
#             @printf("\tNumber of mutations in subclone: %d\n", sresult.output.subclonalmutations[i])
#             @printf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i])
#             @printf("\tTime subclone emerges: %.2f\n", sresult.output.clonetime[i])
#             @printf("\tNumber of divisions: %d\n", sresult.output.Ndivisions[i])
#             @printf("\tAverage number of divisions per cell: %.2f\n", sresult.output.avdivisions[i])
#             @printf("\tPopulation size when subclone emerges: %d\n", sresult.output.cloneN[i])
#             @printf("\tParent of subclone (0 is host): %d\n\n", sresult.output.clonetype[i] - 1)
#         end
#     else
#         @printf("No clones, tumour growth was neutral\n\n")
#     end
# end

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
            @printf("\tEvent time: %.2f\n", 
                sresult.input.siminput.tevent[i])
            @printf("\tAverage frequency: %.2f (std: %.2f)\n", 
                        _get_mean_std(sresult.output, :freq, i)...)
            @printf("\tAverage number of mutations in subclone: %d (std: %.2f)\n", 
                        _get_mean_std(sresult.output, :mutations, i)...)
            @printf("\tAverage population size when subclone emerges: %d (std: %.2f)\n", 
                        _get_mean_std(sresult.output, :N0, i)...)
            @printf("\tAverage time when subclone emerges: %.2f (std: %.2f)\n", 
                        _get_mean_std(sresult.output, :time, i)...)
            @printf("\tDistribution of subclone parents (0 is host): ")
            s = join([@sprintf("%d => %d", parent, freq) 
                        for (parent, freq) in parenttype_dist(sresult.output, i)]
                        , ", ")
            @printf("%s\n\n",s)
        end
    else
        @printf("No clones, tumour growth was neutral\n\n")
    end

end

function _get_mean_std(output, subclonefield, i)
    data = [getfield(subclone, subclonefield)[i]
                for simoutput in output 
                    for subclone in simoutput.subclones]
    return mean(data), std(data)
end

function parenttype_dist(output, i)
    parenttypes = [subclone.parenttype[i] - 1 
                    for simoutput in output 
                        for subclone in simoutput.subclones]
    parenttypes = sort(collect(countmap(parenttypes)))
    return parenttypes
end


# function Base.show(io::IO, sresultlist::Vector{MultiSimulation{T}} where T <: SimulationInput)
#     @printf("%d sets of simulations\n", length(sresultlist))
#     for sresult in sresultlist
#         Base.show(io, sresult)
#     end
# end

function saveinput(filename, sresult::Simulation{BranchingInput})
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
            write(io, @sprintf("Subclone output:\n"))
            for (i,subclone) in enumerate(sresult.output.subclones)
                write(io, @sprintf("Subclone %d: \n", i))
                write(io, @sprintf("Input:"))
                write(io, @sprintf("\tSelection strength, s: %.2f\n", sresult.input.siminput.selection[i]))
                write(io, @sprintf("\tEvent time: %.2f\n", sresult.input.siminput.tevent[i]))
                write(io, @sprintf("Output:\n"))
                write(io, @sprintf("\tFrequency: %.2f\n", subclone.freq))
                write(io, @sprintf("\tNumber of mutations in subclone: %d\n", subclone.mutations))
                write(io, @sprintf("\tFitness advantage: %.2f\n", sresult.input.siminput.selection[i]))
                write(io, @sprintf("\tTime subclone emerges (population doublings): %.2f\n", subclone.time))
                write(io, @sprintf("\tNumber of divisions: %d\n", subclone.Ndivisions))
                write(io, @sprintf("\tAverage number of divisions per cell: %.2f\n", subclone.avdivisions))
                write(io, @sprintf("\tPopulation size when subclone emerges: %d\n", subclone.N0))
                write(io, @sprintf("\tParent of subclone (0 is host): %d\n\n", subclone.parenttype - 1))
            end
        else
            @printf("No clones, tumour growth was neutral\n\n")
        end
        

    end

end

function saveinput(filename, multsim::MultiSimulation, i)
    sresult = get_simulation(multsim, i)
    saveinput(filename, sresult)
end
