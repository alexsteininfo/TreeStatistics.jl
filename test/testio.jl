filename = pwd()*"/testinputsave.json"
@show pwd()
input = MultilevelInput()
saveinput(input, filename)
input2 = loadinput(MultilevelInput, filename)
rm(filename)
@test input == input2
