# MyPMF
Package for calculating potential of mean force (PMF) from data of steered molecular dynamics simulation. This code implements the method proposed in

G. Hummer and A. Szabo, Proc. Nat. Acad. Sci. 98, 3658 (2001).

The authors' approach is based on the Jarzynski equality:

C. Jarzynski, Phys. Rev. Lett. 78, 2690 (1997).

Copyright (c) 2021 Ryuichi Okamoto <<ryubee@gmail.com>>  
License: https://opensource.org/licenses/MIT


## Installation
```sh
julia> ]
pkg> add https://github.com/ryubee1980/MyPMF.jl.git
```

## Usage
# pmf_HS
This function calculate PMF from a set of trajectories of pulled coordinate.
```sh
pmf_HS(traj::Array{Float64,3}, ks::Float64, v::Float64, T::Float64 ; L=500 ::Int64, energy_unit="kcal/mol" ::String)
```

traj[:,:,:] is a 3-dimensional Array of the size K x T x 3, where K is the number of sample trajectories and T is the number of time slices. The time slices must be the same for all samples, traj[i,:,1]=traj[j,:,1] for all i and j.

ks is the spring constant of biasing harmonic potential.
v is the (linear) velocity of the biasing potential.
T is the absolute temperature in units of Kelvin (K).
L is the number of output data points.

The energy_unit must be either "kcal/mol" or "kJ/mol".
The units of length (l), time (T) can be anything, but they should consistently be used for all the variables and parameters. For example, if we set l=nm, T=ps, and E=kJ/mol, then the units of velocity and the spring constant must be [v]=nm/ps and k=kJ/mol/nm^2.

[![Build Status](https://github.com/ryubee1980/MyPMF.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ryubee1980/MyPMF.jl/actions/workflows/CI.yml?query=branch%3Amain)
