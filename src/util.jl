"""
    show(sresult::Simulation)

Print out summary of simulation.
"""
function Base.show(io::IO, sresult::Simulation)
    @printf("===================================================================\n")
    _showtype(sresult)
    @printf("Input parameters: \n")
    @printf("\tMutation rate: %.2f\n", sresult.input.μ)
    _showrates(sresult.input)
    @printf("\tNumber of clonal mutations: %d\n", sresult.input.clonalmutations)
    @printf("\tNumber of subclones: %d\n\n", sresult.input.numclones)
    if sresult.input.numclones > 0
        freqs = subclonefreq(sresult.output)[1]
        for (i, subclone) in enumerate(sresult.output.subclones)
            @printf("Subclone %d \n", i)
            @printf("\tFrequency: %.2f\n", freqs[i])
            @printf("\tNumber of mutations in subclone: %d\n", length(subclone.mutations))
            @printf("\tFitness advantage: %.2f\n", sresult.input.selection[i])
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

function Base.show(io::IO, sresult::MultiSimulation)
    @printf("===================================================================\n")
    @printf("Number of simulations = %d\n\n", length(sresult.output))
    @printf("Input parameters: \n")
    _showrates(sresult.input)
    @printf("\tNumber of clonal mutations: %d\n", sresult.input.clonalmutations)
    @printf("\tNumber of subclones: %d\n\n", sresult.input.numclones)  

    if sresult.input.numclones > 0
        for i in 1:sresult.input.numclones
            @printf("Subclone %d \n", i)
            @printf("\tFitness advantage: %.2f\n", sresult.input.selection[i])
            @printf("\tEvent time: %.2f\n", 
                sresult.input.tevent[i])
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

function Base.show(io::IO, population::Population)
    siminput = population.input
    @printf("===================================================================\n")
    @printf("Multilevel branching module simulation\n\n")
    @printf("Input parameters: \n")
    @printf("\tMax age = %.2f\n", siminput.pop_age)
    _showrates(siminput)
    @printf("\tBranching rate = %.4f\n", siminput.branchrate)
    @printf("\tNew module size = %d\n", siminput.branchinitsize)
    @printf("\tMature module size = %d\n", siminput.modulesize)
    @printf("\tNumber of clonal mutations: %d\n", siminput.clonalmutations)
    @printf("\tNumber of subclones: %d\n\n", siminput.numclones)  

    @printf("Population data\n")
    @printf("\tPopulation age = %d\n", age(population))
    @printf("\tNumber of modules = %d\n\n", length(population))

    # @printf("Clonal mutations per module")

    if siminput.numclones > 0
        for i in 1:siminput.numclones
            @printf("Subclone %d \n", i)
            @printf("\tFitness advantage: %.2f\n", siminput.selection[i])
            @printf("\tEvent time: %.2f\n", 
                siminput.tevent[i])
            @printf("\tAverage frequency: %.2f (std: %.2f)\n", 
                        _get_mean_std(population.output, :freq, i)...)
            @printf("\tAverage number of mutations in subclone: %d (std: %.2f)\n", 
                        _get_mean_std(population.output, :mutations, i)...)
            @printf("\tAverage population size when subclone emerges: %d (std: %.2f)\n", 
                        _get_mean_std(population.output, :N0, i)...)
            @printf("\tAverage time when subclone emerges: %.2f (std: %.2f)\n", 
                        _get_mean_std(population.output, :time, i)...)
            @printf("\tDistribution of subclone parents (0 is host): ")
            s = join([@sprintf("%d => %d", parent, freq) 
                        for (parent, freq) in parenttype_dist(population.output, i)]
                        , ", ")
            @printf("%s\n\n",s)
        end
    else
        @printf("No clones, tumour growth was neutral\n\n")
    end

end

_showtype(sresult::Simulation{BranchingInput}) = @printf("Branching process \n")
_showtype(sresult::Simulation{MoranInput}) = @printf("Moran process \n")
_showtype(sresult::Simulation{BranchingMoranInput}) = @printf("Branching => Moran process\n")
_showtype(sresult::Simulation{MultilevelInput}) = @printf("Branching => Moran proces (multilevel)\n")


function _showrates(siminput::BranchingInput)
    b, d, μ = siminput.b, siminput.d, siminput.μ
    @printf("\tBirth rate of host population: %.2f\n", b)
    @printf("\tDeath rate of host population: %.2f\n", d)
    @printf("\tEffective mutation rate (μ/β): %.2f\n", μ / ((b - d) / b))
end

function _showrates(siminput::MoranInput)
    @printf("\tBirth and death rate of host population: %.2f\n", siminput.bdrate)
end

function _showrates(siminput::Union{BranchingMoranInput, MultilevelInput})
    b, d, μ, bdrate = siminput.b, siminput.d, siminput.μ, siminput.bdrate
    @printf("\tBirth rate of host population (branching): %.2f\n", b)
    @printf("\tDeath rate of host population (branching): %.2f\n", d)
    @printf("\tBirth and death rate of host population (moran): %.2f\n", bdrate)
    @printf("\tMutation rate (μ): %.2f\n", μ)
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
