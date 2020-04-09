# beresheet

Trajectory optimization algorithm for following the trajectory used by the [Beresheet Mission](http://www.visit.spaceil.com/)

## Description

This directory contains several julia files and scripts that use methods from  [TrajectoryOptimization.jl](https://github.com/RoboticExplorationLab/TrajectoryOptimization.jl) to find optimal trajectories in the Earth/Moon system.

## Requirements

* [Julia 1.3+](https://docs.julialang.org/en/v1.3/)

## Installation

Clone the containing repository:

    git clone https://github.com/openlunar/trajectory.git

Download the latest version of julia (version should be at least 1.3) from the [julia website](https://julialang.org/). Instructions for Windows, Mac, or Linux are found on that website.

Once julia is downloaded, you can run it from the terminal with the command `julia`. Make sure the version displayed is greater than 1.3

The `Project.toml` and `Manifest.toml` files in the `trajectory/` repository contain information on all the packages and dependencies used in the julia code. Use julia's [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) `Pkg` to download and precompile them. Make sure you are in the `trajectory/` directory before running the following commands.

```julia
julia> using Pkg
julia> Pkg.activate(".")        # activates the Project.toml file in trajectory/
julia> Pkg.precompile()         # precompiles all packages, may take a bit
julia> Pkg.status()             # displays all packages in Project.toml file
```
Or you can use the built in package manager option:

```julia
julia> ]                        # this sends you to the package manager
(v1.3) pkg> activate .          # activates the Project.toml file in the current directory
(trajectory) pkg> precompile    # note that the header has changed to (trajectory)
(trajectory) pkg> status        # displays all packages in Project.toml file
```
You can exit the built in package manager by pressing `<Backspace>`

Once all your packages are downloaded and precompiled, you're good to go!

*Note: the packages need to be precompiled every time you enter a new julia session*

### Atom

While not essential, I find it convenient to run my julia code from [Atom](https://atom.io/) using the [Juno](http://docs.junolab.org/v0.6/man/installation.html) IDE. It provides a *MATLAB*-like interface with a repl (command line), workspace, debugger, and other features. You can run one line at a time of a julia script using `<Shift>+<Enter>`.

## Structure

+ `trajectory/beresheet/`
  + `BeresheetMoonSunEphem.mat` - mat file containing position and velocity data in ECI for Beresheet, the Moon, and the Sun from **February 22, 2019, 2:18:40 to April 12, 2019, 5:18:00 UTC**, the launch and landing dates of the Beresheet Mission
  + `BeresheetOrbitalElements.jl` - julia script calculating the osculating orbital elements of the Beresheet trajectory with respect to Earth
  + `fullBeresheetTraj.bsp` - binary spk file containing ephemeris data for an **optimized** Beresheet trajectory over the dates given above.
  + `fullBeresheetTraj.mat` - mat file containing the same data as fullBeresheetTraj.bsp
  + `fullBeresheetTraj_nextmonth.bsp` - binary spk file containing ephemeris data for an **optimized** Beresheet trajectory over a timespan one month later than the actual launch.
  + `KKT_Beresheet.jl` - **deprecated** julia script setting up a KKT optimization problem
  + `mission.bsp` - binary spk file containing ephemeris data for a cislunar trajectory created using `pathed_conic.py`.
  + `MissionTrajectory.jl` - julia script used to create new trajectories mimicking Beresheet's at different dates
  + `PiecewiseTrajectory.jl` - julia script used to piece togetether discontinuous trajectories
  + `setupALTRO.jl` - **deprecated** julia script setting up an optimization problem to be fed into the ALTRO solver
  + `TrajectoryPropagation.jl` - julia script with examples of how to propagate initial conditions in different gravity models.
+ `trajectory/utils/`
  + `OptimizationSetup.jl` - julia file that uses commands from [TrajectoryOptimization.jl](https://github.com/RoboticExplorationLab/TrajectoryOptimization.jl) to setup an ILQR optimization problem inside the `optim()` function.
  + `OrbitalDynamics.jl` - julia file that defines several useful dynamics functions, loads [SPICE](https://github.com/JuliaAstro/SPICE.jl) kernels, creates structs containing data on the Earth, the Moon, and the Sun.
+ `trajectory/kernels/` - contains SPICE kernels

## Usage
The commands shown here can be found in the `MissionTrajectory.jl` script.

```julia
include("utils/OrbitDynamics.jl")               # includes useful unction definitions and structs
include("utils/OptimizationSetup.jl"))          # includes optimization function optim()

furnsh(string(pwd(),"/beresheet/fullBeresheetTraj.bsp"))    # load Beresheet kernel
et0b = utc2et("2019-02-22T02:20:48.183")        # actual Beresheet launch date
etfb = utc2et("2019-04-12T05:18:05.187")        # actual Beresheet landing date
boddef("beresheet",-5441)                       # define ephemeris object -5441 as Beresheet
xx_beresheet = [spkezr("beresheet",et,"J2000","none","earth")[1] for et = et0b:60:etfb] # load in Beresheet state data
xx_beresheet = convert2array(xx_beresheet)      # array needs to be fed into optim()

# Change the launch date to a month later
month = 27.321661*86400                         # one month in [sec]
xx,uu = optim(xx_beresheet,et0b+month,etfb+month)   # just one month in the future

# Plot the results
xx_array = convert2array(xx)                    # Easier to plot if in array form
plot(xx_array[1,:],xx_array[2,:],title="Trajectory - 1 Month Later",label="")
```
`xx` is an array of [static arrays](https://github.com/JuliaArrays/StaticArrays.jl) containing the six dimensional state data (position and velocity) of the trajectory

```julia
julia> typeof(xx)
Array{SArray{Tuple{6},Float64,1,6},1}
```
`uu` is an array of [static arrays](https://github.com/JuliaArrays/StaticArrays.jl) containing the three dimensional control (acceleration) necessary to make the trajectory dynamically feasible.

```julia
julia> typeof(uu)
Array{SArray{Tuple{3},Float64,1,3},1}
```
### Changing optimization weights
The Beresheet trajectory has been used as the initial guess for most work with this repository until now. However, the user can use any initial trajectory `xx_init`, initial and final ephemeris times `et0, etf`, and can even change the weights used in the optimization problem to allow for more flexibility. Shown below is an example setting the maginitude of the Q and R matrices to their default values used in the objective function. 

```math
\begin{aligned}
  \min_{x_{0:N},u_{0:N-1}} \quad & \ell_f(x_N) + \sum_{k=0}^{N-1} \ell_k(x_k, u_k, dt) \\
  \textrm{s.t.}            \quad & x_{k+1} = f(x_k, u_k), \\
                                 & g_k(x_k,u_k) \leq 0, \\
                                 & h_k(x_k,u_k) = 0.
\end{aligned}
```

```julia
xx,uu = optim(xx_init, et0, etf, Qmag=1e-2, Rmag = 1e2) 
```

If the optimizer is not converging on a solution, the user can either decrease `Qmag` which will allow the output trajectory to deviate more from the initial guess. Or the user can decrease `Rmag` which will allow a larger control policy (maneuvers) to achieve dynamic feasability. 

## Developers

If you find a bug or wish to make a contribution, use the project's
[github issue tracker](https://github.com/openlunar/trajectory/issues).

## License

Copyright (c) 2019--2020, Open Lunar Foundation.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
