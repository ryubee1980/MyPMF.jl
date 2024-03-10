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
pmf_HS(traj::Array{Float64,3}, ks::Float64, v::Float64, T::Float64 ; L=500 ::Int64, energy_unit="kcal/mol" ::String, J_est=0)
```

traj[:,:,:] is a 3-dimensional Array of the size K x Ts x 3, where K is the number of sample trajectories and Ts is the number of time slices. trj[i,:,1], trj[i,:,2], and trj[i,:,3] are the time, displacement of the steered degree of freedom, and the accumulated work, respectively, of the sample trajectory i. The time and the displacement must start from 0. The time slices must be the same for all samples, traj[i,:,1]=traj[j,:,1] for all i and j.

ks is the spring constant of biasing harmonic potential.
v is the (linear) velocity of the biasing potential.
T is the absolute temperature in units of Kelvin (K).
L is the number of output data points (default value is L=500).

The energy_unit must be either "kcal/mol" or "kJ/mol" (kcal/mol by default).
The units of length (l), time (t) can be anything, but they should consistently be used for all the variables and parameters. For example, if we set l=nm, t=ps, and E=kJ/mol, then the units of velocity and the spring constant must be [v]=nm/ps and [ks]=kJ/mol/nm^2.

If J_est=1, it will also calculate the Jarzynski estimation (default value is 0). That is, $e^{-F/kT}$ is simply estimated as the arithmetic mean of $e^{-w/kT}$. Note that the free energy F estimated by this scheme is equal to the PMF only if the spring constant is large enough (stiff-spring limit).

When big_float=true, the PMF is calculated using BigFloat type. This may be necessary when the samples of accumulated works vary in a wide range, which leads to extremely large or small values of $e^{-w/kT}$.

# pmf_HS_norm
Essentially the same as pmf_HS, but the only difference is that pmf_HS_norm gives the PMF divided by the thermal energy k_BT.

[![Build Status](https://github.com/ryubee1980/MyPMF.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ryubee1980/MyPMF.jl/actions/workflows/CI.yml?query=branch%3Amain)
