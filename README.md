# EPBP

This is an implementation of the EPBP algorithm (see the [*pre-print*](http://arxiv.org/abs/1506.05934)).

This code is released under the MIT license (i.e., feel free to do whatever you want with it) however, please do reach out and let me know if you end up using it as it would be helpful to have feedback on it.

## Disclaimer

This project is still under development [as of June 2015] and should be considered as such. The code is released mainly for demonstration purpose and is not claimed to be optimal but we would be grateful for your feedback.

The release [1DN](https://github.com/tlienart/EPBP) reproduces our original results and should allow one to try the method on an arbitrary graphical model with 1D continuous random variables.


## Contact

- lienart <> stats.ox.ac.uk

## How to use
(section under construction...)

The relevant files are:

- `demo.jl` gives the idea of the "main" file where you can define an MRF, define the potentials and define what you'd like to run (see below)
- `lib_*` define the different functions of interest
-- `_support` a few general purpose functions that are used everywhere
-- `_epbp` the file you're likely to be interested where the EPBP algorithm proper is implemented,
-- `_pbp` implementation of Particle BP (following implementation kindly communicated by A. Ihler & D. Frank)
-- `_ep` implementation of "pure EP"
-- `_lbpd` implementation of loopy BP on a discretized space (taken as "ground truth")
-- `lib_doSim` skeleton for running the different simulations multiple times, with different parameters, etc.

### Definition of an MRF

Two standard Graphical Models can be generated directly: the grid and the chain (see `demo.jl`). For another arbitrary structure, the user can enter a list corresponding to the edges (cf. `demoTree` in `demo.jl`). In practice, the following variables need to be defined:

- `nnodes`: number of nodes
- `nedges`: number of edges
- `edge_list`: a `nedges x 2` matrix where a line corresponds to an edge.

A naive example is the dumbbell graph (1--2) for which we'd have
```
nnodes = 2
nedges = 1
edge_list = [1 2; 2 1]
```
note that the reverse edge need to be explicitly mentioned as well.
