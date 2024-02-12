typedict(x) = Dict{Symbol, Any}(
    fn=>fieldtosave(getfield(x, fn)) for fn in fieldnames(typeof(x))
)

inputdict(x) = push!(typedict(x), :type => string(typeof(x)))

fieldtosave(x) = x
fieldtosave(quiescence::AbstractQuiescence) = inputdict(quiescence)


function saveinput(input, filename)
    filename = filename[end-4:end] == ".json" ? filename : filename * ".json"
    inputdict = typedict(input)
    open(filename, "w") do io
        JSON.print(io, inputdict, 4)
    end
end
