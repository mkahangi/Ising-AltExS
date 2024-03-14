# Ising Model Equation of State (EOS) - C++ Version

## Overview

We present a novel construction of the QCD equation of state (EOS) at finite baryon density. Our work combines a recently proposed resummation scheme for lattice QCD results [[1](https://arxiv.org/abs/2102.06660), [2](https://arxiv.org/abs/2202.05574)] with the universal critical behavior at the QCD critical point [[3](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.034901)]. This allows us to obtain a family of equations of state covering the range $0 \leq \mu_B \leq 700$ MeV and $25$ MeV $\leq T \leq 800$ MeV, which match lattice QCD results near $\mu_B=0$ while featuring a critical point in the 3D Ising model universality class.

The position of the critical point can be chosen within the range accessible to beam-energy scan heavy-ion collision experiments. We parameterize the strength of the singularity and the shape of the critical region using a standard parameter set. We impose stability and causality constraints and discuss the available ranges of critical point parameter choices, finding that they extend beyond earlier parametric QCD EOS proposals. Thermodynamic observables including baryon density, pressure, entropy density, energy density, baryon susceptibility, and speed of sound are presented, covering a wide range in the QCD phase diagram relevant for experimental exploration.

## Physics

### Alternative Expansion Scheme

This expansion scheme provides a relationship between $\chi_1$ and $\chi_2$ with the identity below:

$$
\frac{T}{\mu_B}\chi_1^B(T,\mu_B) = \chi_2(T',0)
$$

where the scaling temperature $T'(T,\mu_B)$ is given by:

$$
T'(T,\mu_B) = T\left[1 + \kappa_2^B(T)\left(\frac{\mu_B}{T}\right)^2 +  ... \right]
$$

<img src="docs/src/mapping.png" alt="Mapping">

### Mapping Ising to QCD

The non-universal mapping from $(r,h) \longleftrightarrow (T, \mu_B)$ is done in two steps. To ensure that the transition of the Ising model $(h=0)$ is aligned with the QCD crossover line, this transformation can be simplified by mapping Ising coordinates $(r,h)$ to alternative scheme coordinates $(T',\mu_B)$ as:

$$
\frac{T' - T_0}{T_C T'_{,T}}  =  -w' h \sin\alpha'_{12}\\
\frac{\mu_B^2 - \mu_{BC}^2}{2\mu_{BC}T_C}  =  w' (-r\rho' - h \cos\alpha'_{12})
$$

#### Merging Ising with Lattice

To ensure that the temperature scaling reproduces lattice results at low $(\mu_B/T)$ and then the critical component from the Ising model only contributes at higher orders in $(\mu_B/T)$, we use:

$$
T'(T,\mu_B) = \underbrace{T'_{\text{lat}}(T,\mu_B)}_{\text{lowest orders in $(\mu_B/T)$ }} + \underbrace{T_{\text{crit}}(T,\mu_B) - \text{Taylor}[T_{\text{crit}}(T,\mu_B)]}_{\text{higher order in $(\mu_B/T)$ }}
$$

where

$$
T'_{\text{crit}}(T,\mu_B) \approx T_0 + \left(\frac{\partial \chi^B_{2,\text{lat}}(T,0)}{\partial T}\Big|_{T_0}\right)^{-1} \frac{n^{\text{crit}}_B(T,\mu_B)}{T^3(\mu_B/T)}f(T,\mu_B)
$$

From the baryon density, we construct all thermodynamic observables using the following:

$$
\begin{align*}
    \frac{\chi_2(T,\mu_B)}{T^2} &= T\left(\frac{\partial ({\color{red}n_B}/T^3)}{\partial \mu_B}\right)\bigg|_{T}\\
    \frac{P(T,\mu_B)}{T^4} &= \chi_{0,\text{lat}}^B(T,0) + \frac{1}{T}\int_0^{\mu_B}d\mu_{B}^\prime{\color{red}n_B(T,\mu_{B}^\prime)}  {T^3} \\
    \frac{s(T,\mu_B)}{T^3} &= \frac{1}{T^3}\frac{\partial P(T,\mu_B)}{\partial T}\Big|_{\mu_B}\\
   \frac{\epsilon(T,\mu_B)}{T^4} &= - \frac{P(T,\hat{\mu}_B)}{T^4} + \frac{s(T,\hat{\mu}_B)}{T^3}  + \hat{\mu}_B\frac{n_B(T,\hat{\mu}_B)}{T^3} \\
    \chi_2^B(T,\mu_B) &= \frac{\partial ({n_B(T,\mu_B)}/T^3)}{\partial \mu_B/T}\bigg|_{T}
\end{align*}
$$

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

### input/parameterfile

Specify the parameter choice:

- `lowT_out`: Minimum Temperature. Example: `lowT_out = 10`
- `highT_out`: Maximum Temperature. Example: `highT_out = 10`
- `T_step`: Stepsize in Temperature. Example: `T_step = 1`
- `lowMU_out`: Minimum baryon chemical potential. Example: `lowMU_out = 0`
- `highMU_out`: Maximum baryon chemical potential. Example: `highMU_out = 700`
- `muB_step`: Stepsize in Baryon chemical potential. Example: `muB_step = 1`
- `muBC`: Location of the critical point in chemical potential
- `alpha12`: Angle between the axis. Example:
