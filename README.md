1`# Ising Model EoS ---C++ VERSION
##  Overview

We present a novel construction of the QCD equation of state (EOS) at finite baryon density. Our work combines a recently proposed resummation scheme  [[2](https://arxiv.org/abs/2102.06660),[3](https://arxiv.org/abs/2202.05574)] for lattice QCD results with the universal critical behavior [[1]](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.034901) at the QCD critical point. This allows us to obtain a family of equations of state in the range $0 \leq \mu_B \leq 700$ MeV and 25 MeV $\leq T \leq 800$ MeV, which match lattice QCD results near $\mu_B=0$ while featuring a critical point in the 3D Ising model universality class.
The position of the critical point can be chosen within the range accessible to beam-energy scan heavy-ion collision experiments. The strength of the singularity and the shape of the critical region are parameterized using a standard parameter set. 
We impose stability and causality constraints and discuss the available ranges of critical point parameter choices, finding that they extend beyond earlier parametric QCD EOS proposals. We present thermodynamic observables, including baryon density, pressure, entropy density, energy density, baryon susceptibility and speed of sound, 
that cover a wide range in the QCD phase diagram relevant for experimental exploration. 

# 1. How this Code Works

## Quick Start

To quickly get started with the Ising Model EOS code, follow these steps:


1. cloning the [Ising repository](https://github.com/mkahangi/Ising-AltExS.git),
  ```git clone https://github.com/mkahangi/Ising-AltExS.git```

2. Navigate to the cloned directory and clean previous builds:
  ```cd Ising-AltExS```
    ``` make clean```
3. Build and run the code using:
  ```make run```
 
After successfully building and running the code, you can customize parameters from the input/parameterfile.
