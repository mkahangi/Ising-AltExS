# Ising Model Equation of State (EOS) - C++ Version

## Overview

We present a novel construction of the QCD equation of state (EOS) at finite baryon density. Our work combines a recently proposed resummation scheme for lattice QCD results [[1][paper1], [2][paper2]] with the universal critical behavior at the QCD critical point [[3][paper3]]. This allows us to obtain a family of equations of state covering the range \(0 \leq \mu_B \leq 700\) MeV and \(25\) MeV \(\leq T \leq 800\) MeV, which match lattice QCD results near \(\mu_B=0\) while featuring a critical point in the 3D Ising model universality class. This work is based on [[4][paper4]]

[//]: # "Link References"
[paper1]: https://arxiv.org/abs/2102.06660
[paper2]: https://arxiv.org/abs/2202.05574
[paper3]: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.034901
[paper4]:  https://arxiv.org/pdf/2402.08636.pdf

...

## How this Code Works

### Quick Start

To quickly get started with the Ising Model EOS code, follow these steps:

1. Clone the [Ising repository](https://gitlab.com/nsf-muses/module-ising-eos/ising_eos):

    ```bash
    git clone git@gitlab.com:nsf-muses/module-ising-eos/ising_eos.git
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

### `input/parameterfile`

Specify the parameter choice:

- `lowT_out`: Minimum Temperature. Example: `lowT_out = 10`
- `highT_out`: Maximum Temperature. Example: `highT_out = 10`
- `T_step`: Stepsize in Temperature. Example: `T_step = 1`
- `lowMU_out`: Minimum baryon chemical potential. Example: `lowMU_out = 0`
- `highMU_out`: Maximum baryon chemical potential. Example: `highMU_out = 700`
- `muB_step`: Stepsize in Baryon chemical potential. Example: `muB_step = 1`
- `muBC`: Location of the critical point in chemical potential
- `alpha12`: Angle between the axis. Example:
