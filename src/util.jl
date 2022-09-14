const multilevel_simulation = runsimulation
const run1simulation = runsimulation
const multilevel_simulation_timeseries = runsimulation_timeseries


function parenttype_dist(output, i)
    parenttypes = [subclone.parenttype[i] - 1 
                    for simoutput in output 
                        for subclone in simoutput.subclones]
    parenttypes = sort(collect(countmap(parenttypes)))
    return parenttypes
end

