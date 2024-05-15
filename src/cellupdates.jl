
"""
    addmutations!(cell::Cell[, cell2::Cell], μ, mutID, rng, mutationdist, Δt=Δt;
        timedepmutationsonly=false)

Add new mutations to `cell` and `cell2` (if it is given) for each mutational process listed
in `mutationdist` (e.g. :poisson, :fixed, :poissontimedep etc) with mean number of mutations
(per division or per unit time) for each process given in vector `μ`.

Time dependent mutations are assigned to both cells, as if the parent had accumulated them
throughout its lifetime (duration `Δt`). Non-time dependent mutations are assigned to each
cell individually. If `timedepmutationsonly` then non-time dependent mutations are not
assigned.
"""
function addmutations! end

function addmutations!(
    cell1::Cell,
    cell2::Cell,
    μ, mutID,
    rng, mutationdist,
    Δt=Δt;
    timedepmutationsonly=false
)
    for (μ0, mutationdist0) in zip(μ, mutationdist)
        if mutationdist0 == :poissontimedep || mutationdist0 == :fixedtimedep
        #for mutations that are time dependent we add the mutations that were accumulated by
        #the parent during its lifetime to the child cells at division
            numbermutations = numbernewmutations(rng, mutationdist0, μ0, Δt=Δt)
            mutID = addnewmutations!(cell1, cell2, numbermutations, mutID)
        elseif !timedepmutationsonly
        #non-time dependent mutations are added to child cells at division
            numbermutations = numbernewmutations(rng, mutationdist0, μ0)
            mutID = addnewmutations!(cell1, numbermutations, mutID)
            numbermutations = numbernewmutations(rng, mutationdist0, μ0)
            mutID = addnewmutations!(cell2, numbermutations, mutID)
        end
    end
    return mutID
end

function addmutations!(
    cell::Cell,
    μ,
    mutID,
    rng,
    mutationdist,
    Δt=Δt;
    timedepmutationsonly=false
)
    for (μ0, mutationdist0) in zip(μ, mutationdist)
        if mutationdist0 == :poissontimedep || mutationdist0 == :fixedtimedep
        #for mutations that are time dependent we add the mutations that were accumulated by
        #the parent during its lifetime to the child cells at division
            numbermutations = numbernewmutations(rng, mutationdist0, μ0, Δt=Δt)
            mutID = addnewmutations!(cell1, numbermutations, mutID)
        elseif !timedepmutationsonly
            #non-time dependent mutations are added to child cells at division
            numbermutations = numbernewmutations(rng, mutationdist0, μ0)
            mutID = addnewmutations!(cell, numbermutations, mutID)
        end
    end
    return mutID
end

"""
    numbernewmutations(rng, mutationdist, μ; Δt=nothing)

Generate the number of new mutations for a given mutational process `mutationdist`
    with mean `μ` (non-time independent, e.g. :fixed, :poisson) or `μΔt` (time dependent,
    e.g. :poissontimedep).
"""
function numbernewmutations(rng, mutationdist, μ; Δt=nothing)
    if mutationdist == :fixed
        return round(Int64, μ)
    elseif mutationdist == :poisson
        return rand(rng, Poisson(μ))
    elseif mutationdist == :geometric
        return rand(rng, Geometric(1/(1+μ)))
    elseif mutationdist == :poissontimedep
        return rand(rng, Poisson(μ*Δt))
    elseif mutationdist == :fixedtimedep
        return round(Int64, μ*Δt)
    else
        error("$mutationdist is not a valid mutation rule")
    end
end

"""
    addnewmutations!(cell::Cell[, cell2::Cell], numbermutations, mutID)

Add `numbermutations` new mutations to `cell` (and identical mutations to `cell2` if given)
starting with id `mutID`. Return `mutID + numbermutations` which will be the next unique
mutation id.
"""
function addnewmutations! end

function addnewmutations!(cell::Cell, numbermutations, mutID)
    #function to add new mutations to cells
    newmutations = mutID:mutID + numbermutations - 1
    append!(cell.mutations, newmutations)
    mutID = mutID + numbermutations
    return mutID
end

function addnewmutations!(cell1::Cell, cell2::Cell, numbermutations, mutID)
    newmutations = mutID:mutID + numbermutations - 1
    append!(cell1.mutations, newmutations)
    append!(cell2.mutations, newmutations)
    mutID = mutID + numbermutations
    return mutID
end

"""
    celldivision!(module, subclones, parentcellid, t, nextID, μ, mutationdist, rng;
        nchildcells=2)

Cell at `module.cells[parentcellid]` divides. If `nchildcells == 1` it is replaced by a
single child cell. If `nchildcells == 2` a second child cell is appended to the end of
`module.cells`. Cell frequencies are update in `subclones`.

Add new mutations according to each mutational process listed in `mutationdist`.
Time-dependent mutations, e.g. `:poissontimedep`, are calculated for the dividing cell,
proportional to its lifetime. Non-time dependent mutations, e.g. `:poisson` are assigned
to the child cell(s) at division.
"""
function celldivision! end

function celldivision!(
    treemodule::TreeModule{T, S},
    subclones,
    parentcellid,
    t,
    nextID,
    μ,
    mutationdist,
    rng;
    nchildcells=2,
    timedepmutationsonly=false
) where {T <: AbstractTreeCell, S}
    alivecells = treemodule.cells
    parentcellnode = alivecells[parentcellid] #get parent cell node
    # assign mutations:
    #    mutations that accumulate with time are assigned to parent cell,
    #    mutations that occur at division are assigned to child cells
    childcellmuts = zeros(Int64, nchildcells)
    for (μ0, mutationdist0) in zip(μ, mutationdist)
        if mutationdist0 == :fixedtimedep || mutationdist0 == :poissontimedep
            Δt = t - parentcellnode.data.latestupdatetime
            parentcellnode.data.mutations += numbernewmutations(
                rng,
                mutationdist0,
                μ0,
                Δt=Δt
            )
            parentcellnode.data.latestupdatetime = t
        elseif !timedepmutationsonly
            for i in 1:nchildcells
                childcellmuts[i] += numbernewmutations(rng, mutationdist0, μ0)
            end
        end
    end
    #create new child cells and add them to alivecells list
    childcell1 = T(
            id=nextID,
            birthtime=t,
            mutations=childcellmuts[1],
            clonetype=parentcellnode.data.clonetype
    )
    #one cell replaces the parent in alivecells
    alivecells[parentcellid] = leftchild!(parentcellnode, childcell1)
    if nchildcells == 2
        childcell2 = T(
            id=nextID + 1,
            birthtime=t,
            mutations=childcellmuts[2],
            clonetype=parentcellnode.data.clonetype
        )
        #the other cell is added to the end of alivecells
        push!(alivecells, rightchild!(parentcellnode, childcell2))
        subclones[parentcellnode.data.clonetype].size += 1
    end
    #adjust subclone sizes
    return treemodule, subclones, nextID + nchildcells
end

function celldivision!(
    cellmodule::CellModule,
    subclones,
    parentcellid,
    t,
    mutID,
    μ,
    mutationdist,
    rng;
    nchildcells=2,
    timedepmutationsonly=false
)
    cellmodule, subclones, Δt = add_new_cells!(
        cellmodule,
        subclones,
        parentcellid,
        t,
        nchildcells
    )
    #add new mutations to both new cells
    if sum(μ) > 0.0
        if nchildcells == 2
            mutID = addmutations!(
                cellmodule.cells[parentcellid],
                cellmodule.cells[end],
                μ,
                mutID,
                rng,
                mutationdist,
                Δt;
                timedepmutationsonly
            )
        else
            mutID = addmutations!(
                cellmodule.cells[parentcellid],
                μ,
                mutID,
                rng,
                mutationdist,
                Δt;
                timedepmutationsonly
            )
        end
    end
    updatetime!(cellmodule, t)
    return cellmodule, subclones, mutID
end

function add_new_cells!(
    cellmodule::CellModule,
    subclones,
    parentcellid,
    t,
    nchildcells
)
    Δt = t - cellmodule.cells[parentcellid].latestupdatetime
    cellmodule.cells[parentcellid].birthtime = t
    cellmodule.cells[parentcellid].latestupdatetime = t
    if nchildcells == 2
        push!(cellmodule.cells, deepcopy(cellmodule.cells[parentcellid]))
        cellmodule.cells[end].id = cellmodule.cells[end-1].id + 1
        cellmodule.cells[end].parentid = cellmodule.cells[parentcellid].id
        subclones[cellmodule.cells[parentcellid].clonetype].size += 1
    end
    return cellmodule, subclones, Δt
end

"""
    cellmutation!(cellmodule, subclones, selectioncoefficient, mutatingcell, t)

Cell mutation occurs in `mutatingcell` to produce a new non-neutral subclone with fitness equal to
`f = 1 + selectioncoefficient`. Implemented so that birth, moran and asymmetric rates are increased
by a factor `f` from the wild-type rates, but death rate remains the same.
"""
function cellmutation!(cellmodule, subclones, selectioncoefficient, mutatingcell, t)
    #add new clone
    subcloneid = length(subclones) + 1
    parentid = getclonetype(mutatingcell)
    wildtype_rates = getwildtyperates(subclones)
    birthrate, deathrate, moranrate, asymmetricrate = get_newsubclone_rates(
        wildtype_rates,
        selectioncoefficient
    )
    newsubclone = Subclone(
        subcloneid,
        parentid,
        t,
        1,
        birthrate,
        deathrate,
        moranrate,
        asymmetricrate
    )
    push!(subclones, newsubclone)
    #change clone type of new cell and update clone sizes
    setclonetype!(mutatingcell, subcloneid)
    if parentid != 0
        subclones[parentid].size -= 1
    end
    return cellmodule, subclones
end

"""
    get_newsubclone_rates(wildtype, selectioncoefficient)
Compute new birth, death, moran and asymmetric rates for a new subclone.

All `wildtype` rates are increased by a factor of `(1 + selectioncoefficient)` (except for death rate
which is unchanged).
"""
function get_newsubclone_rates(wildtype, selectioncoefficient)
    return (
        birthrate = wildtype.birthrate, #wildtype.birthrate * (1 + selectioncoefficient),
        deathrate = wildtype.deathrate + selectioncoefficient, #wildtype.deathrate,
        moranrate = wildtype.moranrate + selectioncoefficient, #wildtype.moranrate * (1 + selectioncoefficient),
        asymmetricrate =  wildtype.asymmetricrate * (1 + selectioncoefficient)
    )
end

"""
    celldeath!(module, subclones, deadcellid, [t, μ, mutationdist, rng])

Cell at `module.cells[deadcellid]` dies and is removed. Cell frequencies are update in
`subclones`. If `module::TreeModule{TreeCell, S}` time-dependent mutations are added
to the dying cell.
"""
function celldeath! end

function celldeath!(
    treemodule::TreeModule,
    subclones::Vector{Subclone},
    deadcellid,
    t,
    μ=nothing,
    mutationdist=nothing,
    rng=nothing
)
    alivecells = treemodule.cells
    deadcellclonetype = alivecells[deadcellid].data.clonetype
    #remove references to dead cell
    killcell!(alivecells, deadcellid, t, μ, mutationdist, rng)
    #remove from alivecell vector
    deleteat!(alivecells, deadcellid)
    #adjust subclone sizes
    subclones[deadcellclonetype].size -= 1
    return treemodule, subclones
end

function celldeath!(cellmodule::CellModule, subclones, deadcell::Integer, args...)
    #frequency of cell type decreases
    clonetype = cellmodule.cells[deadcell].clonetype
    subclones[clonetype].size -= 1
    #remove deleted cell
    deleteat!(cellmodule.cells, deadcell)
    return cellmodule
end



"""
    cellremoval!(module, cells)

Remove `cells` from module without killing them.
"""
function cellremoval!(cellmodule, cells)
    deleteat!(cellmodule.cells, sort!(cells))
    return cellmodule.cells
end

"""
    getnextID(population)
    getnextID(cells)

Get the next usable mutation ID (for Cell type) or cell ID (for AbstractTreeCell type).
"""
getnextID(population::SinglelevelPopulation) = getnextID(population.singlemodule.cells)
getnextID(population::Union{Population, PopulationWithQuiescence}) =
    maximum(map(x->getnextID(x.cells), population))

function getnextID(
    population::Union{Population{T}, PopulationWithQuiescence{T}}
) where T <: TreeModule
    nextID = 1
    for treemodule in population
        for cellnode in treemodule.cells
            if id(cellnode) + 1 > nextID
                nextID = id(cellnode) + 1
            end
        end
    end
    return nextID
end

function getnextID(cells::CellVector)
    if all(no_mutations.(cells))
        return 1
    else
        allmutations = reduce(vcat, [cell.mutations for cell in cells])
        return maximum(allmutations)
    end
end

function getnextID(cells::AbstractTreeCellVector)
    nextID = 1
    for cellnode in cells
        if isnothing(cellnode) continue end
        if id(cellnode) + 1 > nextID
            nextID = id(cellnode) + 1
        end
    end
    return nextID
end

getnextID(cells::TreeCellVector) = maximum(cellnode.data.id for cellnode in cells)
