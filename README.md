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

Simulations are run using the `runsimulation`, `runsimulation_timeseries` or 
`runsimulation_timeseries_returnfinalpop` functions.

To run a simulation an input must be created. Different input types correspond to different 
simulations (see `?SimulationInput` for options). Keyword arguments for each input type and 
default values can be found in the docs for each input type, e.g. see `?BranchingInput` or `MultilevelBranchingMoranInput`. 

By default selection is neutral (i.e. all cells have the 
same fitness). Other selection regimes can be specified by passing an explicit 
`selection::AbstractSelection` argument to `runsimulation`. This can be of type 
`NeutralSelection` (default), `SelectionPredefined` or `SelectionDistribution`.

## Examples 

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

  input = MultilevelBranchingMoranInput(maxmodules=100)
  rng = Random.seed!(12)
  simulation = runsimulation(input, rng)
  ```

- To run a multilevel simulation as above, with new fit mutants occuring at cell division 
  with probability 0.1 with selection coefficients drawn from an Exponential distribution
  with mean 0.2 (maximum of 50 fit mutants):

  ```
  using SomaticEvolution, Random

  input = MultilevelBranchingMoranInput(maxmodules=100)
  selection = SelectionDistribution(Exponential(0.2), 0.2, 50)
  rng = Random.seed!(12)
  simulation = runsimulation(input, selection, rng)
  ```

## Cell implementations

The default implementation relies on storing lists of mutations for each cell, with each
cell having a unique id. A (generally) faster implementation which uses tree structured 
cells is also available. To choose which implementation is used, the cell type can be
passed to runsimulation as the first argument.

- `runsimulation(Cell, input, rng)`: default
- `runsimulation(SimpleTreeCell, input, rng)`: tree-structured cells, dead cells are
  removed from tree.
- `runsimulation(TreeCell, input, rng)`: tree-structured cells, dead cells stay in tree and
  have an additional `alive` field.

## Customisation

Custom simulations can be implemented by defining new `SimulationInput` and/or 
`AbstractSelection` subtypes. This would require defining new `simulate!` and 
`initialize_population` methods that dispatch on the new input/selection types. See source 
code for examples of how these can be implemented.

## Spatial modelling

The second optional argument of `runsimulation` is a type `S<:ModuleStructure` that defaults 
to `WellMixed`. Other types could be implemented to enable spatial arrangement of cells 
within the module.

