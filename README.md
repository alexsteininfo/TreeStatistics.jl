# SomaticEvolution

[![Build Status](https://github.com/jessierenton/SomaticEvolution.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jessierenton/SomaticEvolution.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jessierenton/SomaticEvolution.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jessierenton/SomaticEvolution.jl)

Julia package to simulate single level or multilevel somatic evolution. 

## Installation
To add package to a Julia session run 
```
using Pkg
Pkg.add("https://github.com/jessierenton/SomaticEvolution.jl")
```
## Simulation

To run a simulation first create an input. Different input types correspond to different 
simulations (see `?SimulationInput` for options). Keyword arguments for each input type and default values can be found in the docs for each input type, e.g. see `?BranchingInput`.

- To run a single level branching process simulation, starting from a single cell, until it
  reaches 100 cells: 
  ```   
  using SomaticEvolution, Random

  input = BranchingInput(Nmax=100)
  rng = Random.seed!(12)
  simulation = runsimulation(input, rng)
  ```

- To run a multilevel simulation that starts from a single cell in a single module, grows to 
  a population of 10 modules and then remains in homeostasis:

  ```
  using SomaticEvolution, Random

  input = MultilvelBranchingMoranInput(maxmodules=100)
  rng = Random.seed!(12)
  simulation = runsimulation(input, rng)
  ```
The default implementation relies on storing lists of mutations for each cell, with each
cell having a unique id. A (generally) faster implementation which uses tree structured 
cells is also available. To choose which implementation is used, the cell type can be
passed to runsimulation as the first argument.

- `runsimulation(Cell, input, rng)`: default
- `runsimulation(SimpleTreeCell, input, rng)`: tree-structured cells, dead cells are
  removed from tree.
- `runsimulation(TreeCell, input, rng)`: tree-structured cells, dead cells stay in tree and
  have an additional `alive` field.

Currently, simulations with selection (i.e. `numclones != 0`) are only implemented for 
single-level `Cell` simulations.

## Plots

- To plot the pairwise (fixed) differences between cells or modules use:
  `pairwisedistributionplot(pairwise_fixed_differences(simulation[, idx::Vector{Int}]))`.
  
  If `idx` is given, only those modules are included.
