"""
    multilevel_simulation(IP::InputParameters{MultilevelInput}, rng::AbstractRNG = MersenneTwister; 
        <keyword arguments>)

    Run a simulation process that models somatic evolution in a growing population. Starts 
    with a single cell which grows via a branching process to a "module" of specified size, 
    then switches to a Moran process. New modules are created with given rate, by 
    sampling cells from the parent module.
"""

function multilevel_simulation(IP::InputParameters{MultilevelInput}; 
    rng::AbstractRNG=Random.GLOBAL_RNG)

    populationtracker = initialize_population(
        IP.siminput.Nmax, 
        clonalmutations=IP.siminput.clonalmutations
        )
    population = Population(IP)
    
    for moduletracker in populationtracker
        while true
            moduletracker, newmoduletracker =
                module_simulate_to_branching(
                    moduletracker, 
                    IP,  
                    length(populationtracker) + 1,
                    rng
                )
            if newmoduletracker !== nothing
                push!(populationtracker, newmoduletracker)
            else
                moduletracker, simresults = 
                    processresults!(moduletracker, IP.siminput.Nmax, IP.siminput.numclones, IP.siminput.μ, 
                                IP.siminput.fixedmu, IP.siminput.clonalmutations, IP.ploidy, rng)
            
                #Mimic experimental data by sampling from the true VAF
                sampleddata = 
                    sampledhist(simresults.trueVAF, IP.siminput.Nmax, rng, 
                            detectionlimit = IP.detectionlimit, 
                            read_depth = IP.read_depth, cellularity = IP.cellularity)
                push!(population.output, simresults)
                push!(population.sampled, sampleddata)
                break
            end
        end
    end

    return population

end


function module_simulate_to_branching(moduletracker::ModuleTracker, 
    IP::InputParameters{MultilevelInput}, newmoduleid,
    rng::AbstractRNG=Random.GLOBAL_RNG
    )

    if moduletracker.Nvec[end] < IP.siminput.Nmax
        moduletracker = 
            branchingprocess!(moduletracker, IP.siminput.b, IP.siminput.d, IP.siminput.Nmax, 1, 
                rng, numclones=IP.siminput.numclones, fixedmu=true, selection=IP.siminput.selection, 
                tevent=IP.siminput.tevent, maxclonesize=Inf, tmax=IP.siminput.pop_age)
    end 
    branchtime = 1 / IP.siminput.branchrate .* exptime(rng) + moduletracker.tvec[end]
    if branchtime < IP.siminput.pop_age
        moduletracker = 
            moranprocess!(moduletracker, IP.siminput.bdrate, branchtime, 1, rng::AbstractRNG; 
            numclones = IP.siminput.numclones, fixedmu = true, selection = IP.siminput.selection, 
            tevent = IP.siminput.tevent)
        
        moduletracker, newmoduletracker = 
            sample_new_module!(moduletracker, newmoduleid, IP.siminput.branchinitsize, branchtime, rng)
    
        return moduletracker, newmoduletracker
    else
        return moduletracker, nothing
    end
end

function sample_new_module!(moduletracker, newmoduleid, branchinitsize, branchtime, 
    rng::AbstractRNG)

    sampleids = sample(rng, 1:length(moduletracker.cells), branchinitsize, replace=false)
    newmoduletracker = 
        initializesim_from_cells(moduletracker.cells[sampleids], moduletracker.subclones, 
            newmoduleid, inittime=branchtime)
    moduletracker = celldeath!(moduletracker, sampleids)
    push!(moduletracker.Nvec, moduletracker.Nvec[end] - branchinitsize)
    push!(moduletracker.tvec, branchtime)

    return moduletracker, newmoduletracker
end

function initialize_population(Nmax=nothing; clonalmutations=0)

    #population is defined by ModuleTracker vector with each entry coprresponding to
    #an individual module
    population = ModuleTracker[]

    #initialise population with a single module conisting of a single cell
    initialmodule = initializesim_branching(
        Nmax,
        clonalmutations=clonalmutations
    )
    push!(population, initialmodule)

    return population
end

