filename = pwd()*"/testinputsave.json"
input = MultilevelBranchingInput()
saveinput(input, filename)
input2 = loadinput(MultilevelBranchingInput, filename)
rm(filename)
@test input == input2
