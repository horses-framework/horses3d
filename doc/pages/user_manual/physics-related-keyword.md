title: Physics related keywords
---

[TOC]

## Compressible flow 

| Keyword                   | Description                                                                                                        | Default value          |
|---------------------------|--------------------------------------------------------------------------------------------------------------------|------------------------|
| Mach number               | *REAL*:                                                                                                            | Mandatory keyword     |
| Reynolds number           | *REAL*:                                                                                                            | Mandatory keyword     |
| Prandtl number            | *REAL*:                                                                                                            | 0.72                   |
| Turbulent Prandtl number  | *REAL*:                                                                                                            | Equal to Prandtl       |
| AOA theta                 | *REAL*: Angle of attack (degrees), based on the spherical coordinates polar angle (\(\theta\)) definition           | 0.0                    |
| AOA phi                   | *REAL*: Angle of attack (degrees), based on the spherical coordinates azimuthal angle (\(\varphi\)) definition       | 0.0                    |
| LES model                 | *CHARACTER*(\*): Options are: Vreman, Wale, Smagorinsky, None                                                       | None                   |
| Wall model                | *CHARACTER*:                                                                                                       | linear                 |
| save Lamb vector | *LOGICAL*: Whether to save the Lamb vector or not. The frequency of the saving is determined by `save Lamb vector mode` and `save Lamb vector interval`. Options are: `.true.` or `.false.`. |  `.false.` |
| save Lamb vector mode | *CHARACTER*: The mode to save the Lamb vector. Options are: `time` or `iteration`. |  Mandatory keyword if `save Lamb vector = .true.` |
| save Lamb vector interval | *REAL* or *INTEGER*: The frequency to save the Lamb vector. If `save Lamb vector mode = time`, the Lamb vector is saved every `save Lamb vector interval` seconds of simulated time. If `save Lamb vector mode = iteration`, the Lamb vector is saved every `save Lamb vector interval` iterations. |  Mandatory keyword if `save Lamb vector = .true.` |
| save Lamb vector divergence | *LOGICAL*: Whether to also save the divergence of the Lamb vector, \(\nabla\cdot\mathbf{L}\). It is saved at the same instants as the Lamb vector (and only inside `Lamb vector cube`, if given), in a file `solution_file_<iteration>.divLamb.hsol` with one variable per node. It is computed element by element from the Lamb vector, so it can be discontinuous across element interfaces. It is ignored unless `save Lamb vector = .true.`. Options are: `.true.` or `.false.`. |  `.false.` |
| Lamb vector cube | *REAL*: It is possible to save and compute the Lamb vector only in a cubic region. The six coordinates defining the cube: \([x_0, x_1, y_0, y_1, z_0, z_1]\) |
| save Lamb vector stats | *LOGICAL*: Whether to save the Lamb vector statistics or not. Options are: `.true.` or `.false.`. |  `.false.` |
| save sound velocity squared stats | *LOGICAL*: Whether to save the sound velocity squared statistics or not. A new file `solution_file.SoundVelocitySquared.stats.hsol` is created. Options are: `.true.` or `.false.`. |  `.false.` |
| save gradient sound velocity squared stats | *LOGICAL*: Whether to save the gradient of the sound velocity squared or not. A new file `solution_file.GradientSoundVelocitySquared.stats.hsol` is created. Options are: `.true.` or `.false.`. |  `.false.` |


### Shock capturing 

The shock-capturing module helps stabilize cases with discontinuous solutions, and may also improve the results of under-resolved turbulent cases. It is built on top of a *Sensors* module that detects problematic flow regions, classifying them according to the value of the sensor, \( s \), mapped into the interval \( a \in [0,1] \),


$$
    a = \left\{\begin{array}{ll}
        0, & \text{if } s \leq s_0 - \Delta s / 2, \\
        \frac{1}{2}\left[1+\sin\left(\frac{s-s_0}{\Delta s}\right)\right], & \text{if } s_0 - \Delta s / 2 < s < s_0 + \Delta s / 2,  \\
        1, & \text{elsewhere}.
    \end{array}\right.
$$

The values of \(s_0 = (s_1 + s_2)/2\) and\(\Delta s = s_2 - s_1\) depend on the sensor thresholds\(s_1\) and\(s_2\).

At the moment, flow regions where \(a \leq 0\) are considered smooth and no stabilization algorithm can be imposed there. In the central region of the sensor, with \(0 < a < 1\), the methods shown in the next table can be used and even scaled with the sensor value, so that their intensity increases in elements with more instabilities. Finally, the higher part of the sensor range can implement a different method from the table; however, the intensity is set to the maximum this time.

All the methods implemented introduce artificial dissipation into the equations, which can be filtered with an SVV kernel to reduce the negative impact on the accuracy of the solution. Its intensity is controlled with the parameters\(\mu\) (similar to the viscosity of the Navier-Stokes equations) and\(\alpha\) (scaling of the density-regularization term of the Guermond-Popov flux), which can be set as constants or coupled to the value of the sensor or to a Smagorinsky formulation.

| Keyword | Description | Default value |
|---------|-------------|---------------|
| Enable shock-capturing | *LOGICAL*: Switch on/off the shock-capturing stabilization | .FALSE. |
| Shock sensor | *CHARACTER*: Type of sensor to be used to detect discontinuous regions Options are: <ul><li>Zeros: always return 0</li><li>Ones: always return 1</li><li>Integral: integral of the sensor variable inside each element</li><li>Integral with sqrt: square root of the Integral sensor</li><li>Modal: based on the relative weight of the higher order modes</li><li>Truncation error: estimate the truncation error of the approximation</li><li>Aliasing error: estimate the aliasing error of the approximation</li><li>GMM: clustering sensor based on the divergence of the velocity and the gradient of the pressure</li></ul> | Integral |
| Shock first method | *CHARACTER*: Method to be used in the middle region of the sensor (\(a\in[0,1]\)). Options are: <ul><li>None: Do not apply any smoothing</li><li>Non-filtered: Apply the selected viscous flux without SVV filtering</li><li>SVV: Apply an entropy-stable, SVV-filtered viscous flux</li></ul> | None |
| Shock second method | *CHARACTER*: Method to be used in the top-most region of the sensor (\(a=1\)). Options are: <ul><li>None: Do not apply any smoothing</li><li>Non-filtered: Apply the selected viscous flux without SVV filtering</li><li>SVV: Apply an entropy-stable, SVV-filtered viscous flux</li></ul> | None |
| Shock viscous flux 1 | *CHARACTER*: Viscous flux to be applied in the elements where\(a\in[0,1]\). Options are: <ul><li>Physical</li><li>Guermond-Popov (only with entropy variables gradients)</li></ul> | -- |
| Shock viscous flux 2 | *CHARACTER*: Viscous flux to be applied in the elements where\(a=1\). Options are: <ul><li>Physical</li><li>Guermond-Popov (only with entropy variables gradients)</li></ul> | -- |
| Shock update strategy | *CHARACTER*: Method to compute the variable parameter of the specified shock-capturing approach in the middle region of the sensor. Options are: <ul><li>Constant</li><li>Sensor</li><li>Smagorinsky: only for *non-filtered* and *SVV*</li></ul> | Constant |
| Shock mu 1 | *REAL/CHARACTER(*)*: Viscosity parameter\(\mu_1\), or\(C_s\) in the case of LES coupling | 0.0 |
| Shock alpha 1 | *REAL*: Viscosity parameter\(\alpha_1\) | 0.0 |
| Shock mu 2 | *REAL*: Viscosity parameter\(\mu_2\) | \(\mu_1\) |
| Shock alpha 2 | *REAL*: Viscosity parameter\(\alpha_2\) | \(\alpha_1\) |
| Shock alpha/mu | *REAL*: Ratio between\(\alpha\) and \(\mu\). It can be specified instead of \(\alpha\) itself to make it dependent on the corresponding values of\(\mu\), and it is compulsory when using LES coupling | -- |
| SVV filter cutoff | *REAL/CHARACTER(*)*: Cutoff of the filter kernel,\(P\). If "automatic", its value is adjusted automatically | "automatic" |
| SVV filter shape | *CHARACTER(*)*: Options are: <ul><li>Power</li><li>Sharp</li><li>Exponential</li></ul> | Power |
| SVV filter type | *CHARACTER(*)*: Options are: <ul><li>Low-pass</li><li>High-pass</li></ul> | High-pass |
| Sensor variable | *CHARACTER(*)*: Variable used by the sensor to detect shocks. Options are: <ul><li>rho</li><li>rhou</li><li>rhov</li><li>rhow</li><li>u</li><li>v</li><li>w</li><li>p</li><li>rhop</li><li>grad rho</li><li>div v</li></ul> | rhop |
| Sensor lower limit | *REAL*: Lower threshold of the central sensor region,\(s_1\) | **Mandatory keyword** (except GMM) |
| Sensor higher limit | *REAL*: Upper threshold of the central sensor region,\(s_2\)
| Sensor TE min N      | *INTEGER*: Minimum polynomial order of the coarse mesh used for the truncation error estimation        | 1             |
| Sensor TE delta N    | Polynomial order difference between the solution mesh and its coarser representation                  | 1             |
| Sensor TE derivative | *CHARACTER*: Whether the face terms must be considered in the estimation of the truncation error or not. Options are: <ul><li>Non-isolated</li><li>Isolated</li></ul> | Isolated      |
| Sensor number of clusters | *INTEGER*: Maximum number of clusters to use with the GMM sensor                                       | 2             |
| Sensor min. steps    | *INTEGER*: Minimum number of time steps that an element will remain detected. The last positive value will be used if the sensor "undetects" an element too early | 1             |


#### Spectral vanishing viscosity

The introduction of an SVV-filtered artificial flux helps dissipate high-frequency oscillations. The baseline viscous flux can be chosen as the Navier-Stokes viscous flux or the flux developed by Guermond and Popov. In any case, this flux is expressed in a modal base where it is filtered by any of the following three filter kernels:

- power: \(\hat{F}^{\text{1D}}_i = (i/N)^P\),
- sharp: \(\hat{F}^{\text{1D}}_i = 0\) if \(i \leq P\), \(\hat{F}^{\text{1D}}_i = 1\) elsewhere,
- exponential: \(\hat{F}^{\text{1D}}_i = 0\) if \(i \leq P\), \(\hat{F}^{\text{1D}}_i=\exp\left(-\frac{(i-N)^2}{(i-P)^2}\right)\) elsewhere.

The extension to three dimensions allows the introduction of two types of kernels based on the one-dimensional ones:

- high-pass: \(\hat{F}^{\text{H}}_{ijk} = \hat{F}^{\text{1D}}_i \hat{F}^{\text{1D}}_j \hat{F}^{\text{1D}}_k\),
- low-pass: \(\hat{F}^{\text{L}}_{ijk} = 1 - \left(1-\hat{F}^{\text{1D}}_i\right)\left(1-\hat{F}^{\text{1D}}_j\right)\left(1-\hat{F}^{\text{1D}}_k\right)\),

being the low-pass one more dissipative and, thus, more suited to supersonic cases. The high-pass filter, on the other hand, works better as part of the SVV-LES framework for turbulent cases.

The cutoff parameter \(P\) can be set as "automatic", which uses a sensor to differentiate troubled elements from smooth regions. The stabilisation strategy then depends on the region:

- smooth regions: \(P=4\), \(\mu=\mu_2\), \(\alpha=\alpha_2\),
- shocks: \(P=4\), \(\mu=\mu_1\), \(\alpha=\alpha_1\).

In addition to this, the viscosity \(\mu_1\) can be set to "Smagorinsky" to use the implemented SVV-LES approach. In this case, the \(\mu=\mu_{\text{LES}}\) viscosity is computed following a Smagorinsky formulation with \(C_s=\mu_2\) and the viscosity parameters do not depend on the region anymore,

\[
    \mu = C_s^2 \Delta^2|S|^2, \quad \alpha = \alpha_1.
\]


### Acoustic

The Ffowcs Williams and Hawkings (FWH) acoustic analogy is implemented as a complement to the compressible NS solver. It can run both during the execution of the NS (in-time) or as at a post-process step. The version implemented includes both the solid and permeable surface variations, but both of them for a static body subjected to a constant external flow, i.e. a wind tunnel case scenario. The specifications for the FWH are divided in two parts: the general definitions (including the surfaces) and the observers definitions. The former is detailed in the table below, while the latter are defined in a block section, similar to the monitors (see [monitors](monitors.html)):

```markdown
#define acoustic observer 1
   name     = SomeName
   position = [0.d0, 0.d0, 0.d0]
#end
# end
```

To run the in-time computation, the observers must be defined in the control file. Beware that adding an additional observer will require to run the simulation again. To use the post-process computation, the solution on the surface must be saved at a regular time. Beware that it will need more storage. To run the post-process calculation the horses.tools binary is used, with a control file similar to the one use for the NS simulation (without monitors), and adding the keywords ''tool type'' and ''acoustic files pattern'', as explained in the table below:


| Keyword                        | Description                                                                                                                                                                              | Default value |
|--------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| acoustic analogy               | *CHARACTER(*)*: This is the main keyword for activating the acoustic analogy. The only options is: ''FWH''.                                                                               | --            |
| acoustic analogy permeable     | *LOGICAL*: Defines if uses a permeable or solid approach.                                                                                                                               | .FALSE.       |
| acoustic solid surface         | *CHARACTER(*)*: Array containing the name of each boundary to be used as a surface for integration. In the form: '[bc1,bc2,bc3]'. Mandatory for using the solid surface variant.         | --            |
| acoustic surface file          | *CHARACTER(*)*: Path to a fictitious surface that will be used for integration. It must be tailor-made for the mesh. Mandatory for using the permeable surface variant.             | --            |
| observers flush interval       | *INTEGER*: Iteration interval to flush the observers information to the files.                                                                                                          | 100           |
| acoustic solution save         | *LOGICAL*: Defines whether it saves the NS solution on the surface. Mandatory for post-process computation.                                                                              | .FALSE.       |
| acoustic save timestep        | *REAL*: Controls the time or iteration at which the FWH will be calculated (and saved if specified). If the key is missing it will be done at each timestep.                           | --            |
| acoustic files pattern         | *CHARACTER(*)*: Pattern to the path of all the saved solutions on the surface (To be used in horses.tools for the post-process calculation).                                             | --            |
| tool type                      | *CHARACTER(*)*: Necessary for post-process calculation. Defines the type of post-process of horses.tools. For the FWH analogy the value must be ''fwh post''                          | --            |


## Incompressible Navier-Stokes

Among the various incompressible Navier-Stokes models, `HORSES3D` uses an artificial compressibility 
formulation, which converts the elliptic problem into a hyperbolic system of equations, at the expense of a non divergence–free 
velocity field. However, it allows one to avoid constructing an approximation that satisfies the inf–sup condition. 
This methodology is well suited for use as a fluid flow engine for interface–tracking multiphase flow models, as it allows the 
density to vary spatially. 

The artificial compressibility system of equations is:

\begin{equation}
  \rho_{t} + \vec{\nabla}\cdot\left(\rho\vec{u}\right)=0 ,
  \label{eq:incomressible:continuity}
\end{equation}

\begin{equation}
\left(\rho\vec{u}\right)_t+\vec{\nabla}\cdot\left(\rho \vec{u}\vec{u}\right) = -\vec{\nabla}p + \vec{\nabla}\cdot\left(\frac{1}{\mathrm{Re}}\left(\vec{\nabla}\vec{u} + 
\vec{\nabla}\vec{u}^{T}\right)\right)+\frac{1}{\mathrm{Fr^{2}}}\rho\vec{e}_{g},
\label{eq:incomressible:momentum}
\end{equation}

\begin{equation}
  p_t + \rho_0 c_0^2 \vec{\nabla}\cdot\vec{u} = 0,
  \label{eq:incomressible:ACM}
\end{equation}

The factor \(\rho_0\) is computed as \(max(\rho_1/\rho_1,\rho_2/\rho_1)\).

This solver is run with the binary horses3d.ins. 

| Keyword                          | Description                                                                                       | Default value |
|----------------------------------|---------------------------------------------------------------------------------------------------|---------------|
| reference velocity (m/s)        | *REAL*: Reference value for velocity                                                              | --            |
| number of fluids (1/2)          | *INTEGER*: Number of fluids present in the simulation                                              | 1             |
| maximum density (kg/m^3)        | *REAL*: Maximum value used in the limiter of the density                                           | Huge(1.0)     |
| minimum density (kg/m^3)        | *REAL*: Minimum value used in the limiter of the density                                           | -Huge(1.0)    |
| artificial compressibility factor | *REAL*: Artificial compressibility factor \(c_0^2\)                                                | --            |
| gravity acceleration (m/s^2)    | *REAL*: Value of gravity acceleration                                                             | --            |
| gravity direction                | *REAL*: Array containing direction of gravity. Eg. \( [0.0,-1.0,0.0] \)                             | --            |
| save Lamb vector | *LOGICAL*: Whether to save the Lamb vector or not. The frequency of the saving is determined by `save Lamb vector mode` and `save Lamb vector interval`. Options are: `.true.` or `.false.`. |  `.false.` |
| save Lamb vector mode | *CHARACTER*: The mode to save the Lamb vector. Options are: `time` or `iteration`. |  Mandatory keyword if `save Lamb vector = .true.` |
| save Lamb vector interval | *REAL* or *INTEGER*: The frequency to save the Lamb vector. If `save Lamb vector mode = time`, the Lamb vector is saved every `save Lamb vector interval` seconds of simulated time. If `save Lamb vector mode = iteration`, the Lamb vector is saved every `save Lamb vector interval` iterations. |  Mandatory keyword if `save Lamb vector = .true.` |
| save Lamb vector divergence | *LOGICAL*: Whether to also save the divergence of the Lamb vector, \(\nabla\cdot\mathbf{L}\). It is saved at the same instants as the Lamb vector (and only inside `Lamb vector cube`, if given), in a file `solution_file_<iteration>.divLamb.hsol` with one variable per node. It is computed element by element from the Lamb vector, so it can be discontinuous across element interfaces. It is ignored unless `save Lamb vector = .true.`. Options are: `.true.` or `.false.`. |  `.false.` |
| save Lamb vector stats | *LOGICAL*: Whether to save the Lamb vector statistics or not. Options are: `.true.` or `.false.`. |  `.false.` |



The incompressible Navier Stokes solver has two modes: with 1 fluid and with 2 fluids. The required keywords are listed below.

Table 1: Keywords for incompressible Navier-Stokes solver. Mode with 1 fluid

| Keyword             | Description                        | Default value |
|---------------------|------------------------------------|---------------|
| density (kg/m^3)    | *REAL*: Density of the fluid      | --            |
| viscosity (pa.s)    | *REAL*: Viscosity of the fluid    | --            |

Table 2: Keywords for incompressible Navier-Stokes solver. Mode with 2 fluids

| Keyword                | Description                           | Default value |
|------------------------|---------------------------------------|---------------|
| fluid 1 density (kg/m^3) | *REAL*: Density of the fluid 1       | --            |
| fluid 1 viscosity (pa.s) | *REAL*: Viscosity of the fluid 1     | --            |
| fluid 2 density (kg/m^3) | *REAL*: Density of the fluid 2       | --            |
| fluid 2 viscosity (pa.s) | *REAL*: Viscosity of the fluid 2     | --            |


## Multiphase

The multiphase flow solver implemented in `HORSES3D` is constructed by a combination of the diffuse interface model of Cahn–Hilliard
with the incompressible Navier–Stokes equations with variable density and artificial compressibility. 
This model is entropy stable and guarantees phase conservation with an accurate representation of surface tension effects. 
The modified entropy-stable version approximates:

\begin{equation}
c_t + \vec{\nabla}\cdot\left(c\vec{u}\right) = M_0 \vec{\nabla}^2 \mu,
\label{eq:governing:cahn--hilliard}
\end{equation}

\begin{equation}
\sqrt{\rho}\left(\sqrt{\rho}\vec{u}\right)_t+\vec{\nabla}\cdot\left(\frac{1}{2}\rho \vec{u}\vec{u}\right) 
+\frac{1}{2}\rho\vec{u}\cdot\vec{\nabla}\vec{u}+c\vec{\nabla}\mu
= -\vec{\nabla}p + \vec{\nabla}\cdot\left(\eta\left(\vec{\nabla}\vec{u} + 
\vec{\nabla}\vec{u}^{T}\right)\right)+\rho\vec{g},
\label{eq:governing:momentum-skewsymmetric-sqrtRho}
\end{equation}

\begin{equation}
  p_t + \rho_0 c_0^2 \vec{\nabla}\cdot\vec{u} = 0,
  \label{eq:governing:ACM}
\end{equation}

where \(c\) is the phase field parameter, \(M_0\) is the mobility, \(\mu\) is the chemical potential, 
\(\eta\) is the viscosity and \(c_0\) is the artificial speed of sound. The factor \(\rho_0\) is computed as \(max(\rho_1/\rho_1,\rho_2/\rho_1)\).
Mobility \(M_0\) is computed from the control file parameters chemical characteristic time \(t_{CH}\), interface width \(\epsilon\) and interface tension \(\sigma\) with the formula 
\(M_0 = L_{ref}^2 \epsilon /(t_{CH} \sigma)\). 

The term \(M_0 \vec{\nabla}^2 \mu\) can be implicity integrated to reduce the stiffnes of the problem with the keyword time integration = IMEX. This is only recomended if the value of \(M_0\) is very high so that the time step of the explicit scheme is very small. 

This solver is run with the binary horses3d.mu. 

| Keyword                          | Description                                                                                         | Default value |
|----------------------------------|-----------------------------------------------------------------------------------------------------|---------------|
| fluid 1 density (kg/m^3)        | *REAL*: Density of fluid 1                                                                          | --            |
| fluid 1 viscosity (pa.s)         | *REAL*: Viscosity of fluid 1                                                                        | --            |
| fluid 2 density (kg/m^3)        | *REAL*: Density of fluid 2                                                                          | --            |
| fluid 2 viscosity (pa.s)         | *REAL*: Viscosity of fluid 2                                                                        | --            |
| reference velocity (m/s)         | *REAL*: Reference value for velocity                                                                 | --            |
| maximum density (kg/m^3)        | *REAL*: Maximum value used in the limiter of the density                                            | Huge(1.0)     |
| minimum density (kg/m^3)        | *REAL*: Minimum value used in the limiter of the density                                            | -Huge(1.0)    |
| artificial compressibility factor | *REAL*: Artificial compressibility factor \({c_0}^2\)                                                 | --            |
| gravity acceleration (m/s^2)     | *REAL*: Value of gravity acceleration                                                              | --            |
| gravity direction                 | *REAL*: Array containing direction of gravity. Eg. \([0.0, -1.0, 0.0]\)                               | --            |
| velocity direction               | *REAL*: Array containing direction of velocity used for the outflow BC. Eg. \([1.0, 0.0, 0.0]\)       | --            |
| chemical characteristic time (s) | *REAL*: \(t_{CH}\) controls the speed of the phase separation                                         | --            |
| interface width (m)              | *REAL*: \(\epsilon\) controls the interface width between the phases                                   | --            |
| interface tension (N/m)          | *REAL*: \(\sigma\) controls the interface tension between the phases                                   | --            |
| save Lamb vector | *LOGICAL*: Whether to save the Lamb vector or not. The frequency of the saving is determined by `save Lamb vector mode` and `save Lamb vector interval`. Options are: `.true.` or `.false.`. |  `.false.` |
| save Lamb vector mode | *CHARACTER*: The mode to save the Lamb vector. Options are: `time` or `iteration`. |  Mandatory keyword if `save Lamb vector = .true.` |
| save Lamb vector interval | *REAL* or *INTEGER*: The frequency to save the Lamb vector. If `save Lamb vector mode = time`, the Lamb vector is saved every `save Lamb vector interval` seconds of simulated time. If `save Lamb vector mode = iteration`, the Lamb vector is saved every `save Lamb vector interval` iterations. |  Mandatory keyword if `save Lamb vector = .true.` |
| save Lamb vector divergence | *LOGICAL*: Whether to also save the divergence of the Lamb vector, \(\nabla\cdot\mathbf{L}\). It is saved at the same instants as the Lamb vector (and only inside `Lamb vector cube`, if given), in a file `solution_file_<iteration>.divLamb.hsol` with one variable per node. It is computed element by element from the Lamb vector, so it can be discontinuous across element interfaces. It is ignored unless `save Lamb vector = .true.`. Options are: `.true.` or `.false.`. |  `.false.` |
| save Lamb vector stats | *LOGICAL*: Whether to save the Lamb vector statistics or not. Options are: `.true.` or `.false.`. |  `.false.` |


## Particles

`Horses3d` includes a two-way coupled Lagrangian solver.
Particles are tracked along their trajectories, according to the simplified particle equation of motion, where only contributions from Stokes drag and gravity are retained,
\begin{equation}
\label{eq:part_motion}
\frac{d y_i}{dt} = u_i, \quad \frac{d u_i}{dt} = \frac{v_i - u_i}{\tau_p} + g_i,
\end{equation}
where \(u_i\) and \(y_i\) are the \emph{ith} components of velocity and position of the particle, respectively. Furthermore, \(v_i\) accounts for the continuous velocity of the fluid at the position of the particle.  We consider spherical Stokesian particles, so their mass and aerodynamic response time are \(m_p = \rho_p \pi D_p^3/6\) and \(\tau_p = \rho_p D_p^2 / 18\mu \), respectively, \(\rho_p\) being the particle density and \(D_p\) the particle diameter. 

Each particle is considered to be subject to a radiative heat flux \(I_o\). The carrier phase is transparent to radiation, whereas the incident radiative flux on each particle is completely absorbed. Because we focus on relatively small volume fractions, the fluid-particle medium is considered to be optically thin. Under these hypotheses, the direction of the radiation is inconsequential, and each particle receives the same radiative heat flux, and its temperature \(T_p\) is governed by
\begin{equation}
\label{eq:part_energy}
\frac{d}{dt} (m_p c_{V,p} T_p) = \frac{\pi D_p^2}{4} I_o - \pi D_p^2 h (T_p-T),
\end{equation}
where \(c_{V,p}\) is the specific heat of the particle, which is assumed to be constant with respect to temperature. \(T_p\) is the particle temperature and \(h\) is the convective heat transfer coefficient, which for a Stokesian particle can be calculated from the Nusselt number \(Nu = hD_p/k = 2\).

In practical simulations, integrating the trajectory of every particle is too expensive. Therefore, particles are agglomerated into parcels, each of them accounting for many particles with the same physical properties, position, velocity, and temperature. The evolution of the parcels is tracked with the same set of equations presented for the particles.

The two-way coupling means that fluid flow is modified because of the presence of particles. Therefore, the Navier-Stokes equations are enriched with the following source terms:

\begin{equation}
\boldsymbol{S} = \beta\left[\begin{array}{c} 0 \\
                                                                       \sum_{n=1}^{N_p} \frac{m_p}{\tau_p} (u_{1,n}-v_1)\delta(\mathbf{x} - \mathbf{y}_n)  \\
                                                                       \sum_{n=1}^{N_p} \frac{m_p}{\tau_p} (u_{2,n}-v_2)\delta(\mathbf{x} - \mathbf{y}_n) \\
                                                                      \sum_{n=1}^{N_p} \frac{m_p}{\tau_p} (u_{3,n}-v_3)\delta(\mathbf{x} - \mathbf{y}_n) \\
                                                                      \sum_{n=1}^{N_p} \pi D_p^2 h (T_{p,n} - T) \delta(\mathbf{x} - \mathbf{y}_n )
\end{array}\right],
\label{eq:particlessource}
\end{equation}

where \(\delta\) is the Dirac delta function, \(N_p\) is the number of parcels, \(\beta\) is the number of particles per parcel and \(u_{i,n}\), \(\mathbf{y}_{i,n}\), \(T_{p,n}\) are the velocity, spatial coordinates, and temperature of the parcel \emph{nth}.
The dimensionless form of the Navier Stokes equations can be seen in the appendix at the end of this document. 

Particles are solved in a box domain inside the flow domain. The box is defined with the keywords ``minimum box'' and ``maximum box''. The boundary conditions for the particles are defined with the keyword ``bc box''. Possible options are inflow/outflow, periodic and wall (with perfect rebound).
 
| Keyword                          | Description                                                                                      | Default value |
|----------------------------------|--------------------------------------------------------------------------------------------------|---------------|
| lagrangian particles             | *LOGICAL*: If .true. activates particles                                                         | .false.       |
| stokes number                    | *REAL*: Stokes number which for Stokesian particles is \(St=\frac{\rho_p D_p^2 u_o}{18 L_o \mu_o}\) | --            |
| Gamma                            | *REAL*: Ratio between specific heat of particles and fluid \(\Gamma=c_{v,p}/c_{v}\)                | --            |
| phi_m                            | *REAL*: Ratio between total mass of particles and fluid \(\phi_m=\frac{m_p N_p}{\rho_o L_o^3}\)    | --            |
| Radiation source                 | *REAL*: Non-dimensional radiation source intensity \(I_o^*=\frac{I_o D_p}{4k_oT_o}\)                | --            |
| Froude number                    | *REAL*: Froude number \(Fr=\frac{u_o}{\sqrt{g_o L_o}}\)                                            | --            |
| high order particles source term | *LOGICAL*: Source term with high order Dirac delta or averaged in the whole element               | .false.       |
| number of particles              | *INTEGER*: Total number of parcels in the simulation                                              | --            |
| particles per parcel             | *REAL*: \(\beta\) particles per parcel                                                              | --            |
| Gravity direction                | *INTEGER*: Array with direction of gravity. Only required if Fr number is specified. [0,0,-1]     | --            |
| particles file                   | *CHARACTER(\*)*: Path to file with initial position of the particles.                             | --            |
| vel and temp from file           | *LOGICAL*: If .true. Initial velocity and temperature of particles read from file.               | --            |
| injection                        | *LOGICAL*: If .true. injection of particles through a face of the box.                             | --            |
| particles injection              | *INTEGER*: Array with a vector indicating the direction of the injection. Eg., [0,1,0]            | --            |
| particles per step               | *INTEGER*: Number of particles injected per time step.                                            | --            |
| particles iter period            | *INTEGER*: Iteration period for particles injection. Set to 1 to inject particles every time step.| --            |
| particles injection velocity     | *REAL*: Array with particles injection non-dimensional velocity. Eg., [0.d0,1.d0,0.d0]            | --            |
| particles injection temperature  | *REAL*: Particles injection non-dimensional temperature                                           | --            |
| minimum box                      | *REAL*: Array with minimum x,y,z coordinates of box with particles. Eg., [0.d0,0.d0,0.d0]         | --            |
| maximum box                      | *REAL*: Array with maximum x,y,z coordinates of box with particles. Eg., [4.d-2,1.6d-1,4.d-2]     | --            |
| bc box                           | *INTEGER*: Array with boundary conditions for particles box in the form [yz,xz,xy].               | --            |

@warning 
- Lagrangian particles are only implemented for the compressible Navier Stokes
- Lagrangian particles do not support MPI
@endwarning

## Complementary Modes

### Wall functions

The wall function overwrites the viscous flux on the specified boundaries based on an specific law using a Newman condition. It must be used as a complement of no slip boundary condition. The table below shows the parameters that can be set in the control file. The frictional velocity is calculated using the instantaneous values of the first node (either Gauss or Gauss-Lobatto) of the element neighbour of the face element (at the opposite side of the boundary face). Currently is only supported for the compressible Navier-Stokes solver.

The standard wall function uses the Reichardt law, solving the algebraic non-linear equation using the newton method to obtain the frictional velocity. The ABL function uses the logarithmic atmospheric boundary layer law, using the aerodynamic roughness; the frictional velocity is without using any numerical method.

| Keyword                        | Description                                                                                      | Default value |
|--------------------------------|--------------------------------------------------------------------------------------------------|---------------|
| Wall Function                  | *CHARACTER(\*)*: This is the main keyword for activating the wall function. Identifies the wall law to be used. Options are:
                                  - Standard: uses the Reichardt law.
                                  - ABL: uses the atmospheric boundary layer law.                                      | --            |
| Wall Function Boundaries       | *CHARACTER(\*)*: Array containing the name of each boundary to be used. In the form: '[bc1,bc2,bc3]'. Mandatory for using the wall function.                   | --            |
| Wall Function kappa            | *REAL*: von Karman constant                                                                     | 0.38          |
| Wall Function C                | *REAL*: Log law 'C' constant                                                                    | 4.1           |
| Wall Function Seed             | *REAL*: Initial value for the newton method                                                      | 1.0           |
| Wall Function Damp             | *REAL*: Initial value damp for the newton method                                                | 1.0           |
| Wall Function Tolerance        | *REAL*: Tolerance for the newton method                                                          | \(10^{-10}\)    |
| Wall Function max iter         | *INTEGER*: Maximum number of iterations for the newton method                                    | 100           |
| Wall Roughness                 | *REAL*: Aerodynamic roughness for the ABL wall function. Mandatory value for the ABL law.        | --            |
| Wall Plane Displacement        | *REAL*: Plane displacement due to roughness for the ABL wall function                             | 0.0           |
| Wall Function Use Average      | *LOGICAL*: Use the time average of the velocity in the wall function, each time step the time average is recalculated.                                       | .FALSE.       |

### Tripping 

A numerical source term is added to the momentum equations to replicate the effect of a tripping mechanism used commonly in explerimental tests. The forcing is described via the product of two independent functions: one that depends streamwise and vertical directions (space only) and the other one describing the spanwise direction and time (space and time).
It can be used for the compressible NS, both LES and RANS.
The keywords for the trip options are:

| Keyword              | Description                                                                                                              | Default value |
|----------------------|--------------------------------------------------------------------------------------------------------------------------|---------------|
| use trip             | *LOGICAL*: This is the main keyword for activating the trip                                                              | .FALSE.       |
| trip time scale      | *REAL*: Time interval between the change of the time dependent part of the trip.                                        | **Mandatory** |
| trip number of modes | *INTEGER*: Number of Fourier modes in the spanwise direction of the trip.                                                | **Mandatory** |
| trip z points        | *INTEGER*: Number of points to create the Fourier Transformation of the spanwise direction, it must be greater than the number of modes and should be ideally equal to the number of discretization points of the mesh in the same direction. | **Mandatory** |
| trip attenuation     | *REAL ARRAY(2)*: Length scale of the gaussian attenuation of the trip, the first position is the streamwise direction and the second is the wall-normal direction. | **Mandatory** |
| trip zone            | *CHARACTER(\*) ARRAY(:)*: Boundary condition name that constrains at least one surface where the trip center is located. It can be either one or two boundary conditions, the latter used to generate a trip in two different positions (i.e. pressure and suction sides of an airfoil). | **Mandatory** |
| trip center          | *REAL*: Position of the origin of the trip in the streamwise direction.                                                  | **Mandatory** |
| trip center 2        | *REAL*: Position of the origin of the second trip, if used, in the streamwise direction.                                 | --            |
| trip amplitude       | *REAL*: Maximum time varying amplitude of the trip.                                                                      | 1.0           |
| trip amplitude steady| *REAL*: Maximum steady amplitude of the trip.                                                                             | 0.0           |
| random seed 1        | *INTEGER*: Number used to initialize the random number generator of the trip. It can vary in different simulations but must remain constant for a restart. | 930187532     |
| random seed 2        | *INTEGER*: Number used to initialize the random number generator of the trip. It can vary in different simulations but must remain constant for a restart. | 597734650     |

## Acoustics

| Keyword                   | Description                                                                                                        | Default value          |
|---------------------------|--------------------------------------------------------------------------------------------------------------------|------------------------|
| base solver  | *CHARACTER*: Specifies which solver was run to compute the base flow. Options are: `ns`/`ins`/`mu`. | Mandatory keyword |
| qbase | *CHARACTER*: Whether to set a uniform base flow or read it from a file. Options are: `file` or `uniform`. | Mandatory keyword |
| qbase vector | *REAL*: The vector of the base flow. The order is: mean density, mean velocity, mean pressure, mean sound velocity squared. | Mandatory keyword if `qbase = uniform` |
| qbase file name | *CHARACTER*: The path to the stats file of the base flow. | Mandatory keyword if `qbase = file`. |
| use source term | *LOGICAL*: Whether to compute the source term or not. Options are: `.true.` or `.false.`. | `.false.` |
| APE number | *INTEGER*: The number of the APE to use. Options are: 4. | 4 |
| sound velocity squared base file name | *CHARACTER*: The path to the stats file of the base sound velocity squared. | Mandatory keyword if `base solver = ns`. |
| gradient sound velocity squared base file name | *CHARACTER*: The path to the stats file of the base gradient sound velocity squared. | Mandatory keyword if `base solver = ns`. |
| use Lamb vector | *LOGICAL*: Whether to compute the Lamb vector in the source term or not. Options are: `.true.` or `.false.`. | `.false.` |
| Lamb vector base | *CHARACTER*: Whether to set a uniform base Lamb vector or read it from a file. Options are: `file` or `uniform`. | Mandatory keyword if `use Lamb vector = .true.` |
| Lamb vector base vector | *REAL*: The vector of the base Lamb vector. | Mandatory keyword if `Lamb vector base = uniform` |
| Lamb vector base file name | *CHARACTER*: The path to the stats file of the base Lamb vector. | Mandatory keyword if `Lamb vector base = file`. |
| save pressure squared stats | *LOGICAL*: Whether to save the squared pressure statistics or not. Options are: `.true.` or `.false.`. | `.false.`. |


When the base solver is compressible Navier-Stokes or multiphase, it is possible to input either the Mach number or a reference velocity. In the case of compressible Navier-Stokes, the Mach number and the reference velocity are related by
\begin{equation}
V_{\text{ref}} = \text{M} \; \sqrt{ \gamma R T_{\text{ref}}},
\end{equation}
where \(\gamma\) is the heat ratio, \(R\) is the gas constant, and \(T_{\text{ref}}\) is a reference temperature.
In the case of multiphase solver, the Mach number and the reference velocity are related by
\begin{equation}
V_{\text{ref}} = \text{M} \; a^2_1,
\end{equation}
where \(a_1\) is the sound velocity in the first fluid.
When the base solver is incompressible Navier-Stokes, it is only possible to input a reference velocity.

### Lamb vector time interpolation

The data of the Lamb vector from the base flow simulation (NS, iNS or MU) is stored in different `.Lamb.hsol` files at several time instants. The time steps of the acoustics simulation might not match these time instants. Thus, a time interpolation is needed.

The Lamb vector data is stored in files with names: `caseName_00000.Lamb.hsol`, `caseName_00001.Lamb.hsol`, `caseName_00002.Lamb.hsol`, etc. To perform the time interpolation, we consider the files `dir/to/lamb/caseName*`, where `Lamb vector dir = "dir/to/lamb/"` and `Lamb vector file name = "caseName"`. Next, we read all the times from these files and sort them in ascending order. Then, at each time step of the acoustics simulation, we interpolate in time.

When `Lamb vector time interpolation = "constant"`, the Lamb vector is assumed constant throughout all the simulation and it is equal to the value of the Lamb vector from the first file. When `Lamb vector time interpolation = "linear"`, the Lamb vector is interpolated linearly between the previous and next data.

| Keyword                   | Description                                                                                                        | Default value          |
|---------------------------|--------------------------------------------------------------------------------------------------------------------|------------------------|
| Lamb vector time interpolation | *CHARACTER*: The type of time interpolation of the Lamb vector. Options are: `constant` or `linear`. | `constant` |
| Lamb vector dir | *CHARACTER*: The path to the directory containing the Lamb vector files. | Mandatory keyword if `use Lamb vector = .true.`. |
| Lamb vector file name | *CHARACTER*: The basename of the files containing the Lamb vector data. | Mandatory keyword if `use Lamb vector = .true.`. |

### Lamb vector inside region

It is possible to compute the contribution of the Lamb vector to the source term only in a region of the domain (a cube). The user specifies the six corners defining the straight-sided cube with the keyword `Lamb vector cube`. The default behavior is to compute the Lamb vector in the whole domain.

| Keyword                   | Description                                                                                                        | 
|---------------------------|--------------------------------------------------------------------------------------------------------------------|
| Lamb vector cube | *REAL*: The six coordinates defining the cube: \([x_0, x_1, y_0, y_1, z_0, z_1]\) |