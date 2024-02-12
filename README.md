# Ising Model Equation of State (EOS) - C++ Version

## Overview

This repository presents a novel construction of the Quantum Chromodynamics (QCD) equation of state (EOS) at finite baryon density. Our work combines a recently proposed resummation scheme [1](https://arxiv.org/abs/2102.06660) [2](https://arxiv.org/abs/2202.05574) for lattice QCD results with the universal critical behavior [3](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.034901) at the QCD critical point. This allows us to obtain a family of equations of state in the range \(0 \leq \mu_B \leq 700\) MeV and \(25\) MeV \(\leq T \leq 800\) MeV, matching lattice QCD results near \(\mu_B=0\) while featuring a critical point in the 3D Ising model universality class.

The position of the critical point can be chosen within the range accessible to beam-energy scan heavy-ion collision experiments. We parameterize the strength of the singularity and the shape of the critical region using a standard parameter set. Stability and causality constraints are imposed, and we discuss the available ranges of critical point parameter choices, finding that they extend beyond earlier parametric QCD EOS proposals. The repository provides thermodynamic observables, including baryon density, pressure, entropy density, energy density, baryon susceptibility, and speed of sound, covering a wide range in the QCD phase diagram relevant for experimental exploration.

## 1. How this Code Works

### Quick Start

To quickly get started with the Ising Model EOS code, follow these steps:

1. Clone the [Ising repository](https://github.com/mkahangi/Ising-AltExS.git):
    ```bash
    git clone https://github.com/mkahangi/Ising-AltExS.git
    ```

2. Navigate to the cloned directory and clean previous builds:
    ```bash
    cd Ising-AltExS
    make clean
    ```

3. Build and run the code using:
    ```bash
    make run
    ```

After successfully building and running the code, you can customize parameters from the `input/parameterfile`.

### input/parameterfile

Specify the parameter choice:

- `lowT_out`: Minimum Temperature. Example: `lowT_out = 10`
- `highT_out`: Maximum Temperature. Example: `highT_out = 10`
- `T_step`: Stepsize in Temperature. Example: `T_step = 1`
- `lowMU_out`: Minimum baryon chemical potential. Example: `lowMU_out = 0`
- `highMU_out`: Maximum baryon chemical potential. Example: `highMU_out = 700`
- `muB_step`: Stepsize in Baryon chemical potential. Example: `muB_step = 1`
- `muBC`: Location of the critical point in chemical potential
- `alpha12`: Angle between the axis. Example: `alpha12=90`
- `ww`: Scaling parameter. Example: `ww = 2`
- `rho`: Scaling parameter. Example: `rho = 2`
