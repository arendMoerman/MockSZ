---
title: "Fast generation of mock Sunyaev-Zel'dovich maps"
tags:
  - Sunyaev-Zel'dovich effect
  - Galaxy clusters
  - Python
  - C++
authors:
  - name: Arend Moerman
    orcid: 0000-0002-0475-6134
    corresponding: true
    affiliation: 1
affiliations:
  - name: Faculty of Electrical Engineering, Mathematics and Computer Science, Delft University of Technology, Mekelweg 4, 2628 CD, Delft, The Netherlands
    index: 1
date: 26 April 2023
bibliography: paper.bib
---

# Summary
Galaxy clusters are some of the largest known gravitationally bound structures in the Universe.
Studying these structures 
The plasma in galaxy clusters is known to interact with the cosmic microwave background (CMB), leaving a characteristic spectral distortion in the observed CMB spectrum, which is strongest in the mm/sub-mm. 
This particular 

`MockSZ` is a Python interface for fast calculations of the thermal and kinematic SZ effects from galaxy clusters. 
It calculates the SZ effect by solving the collisionless Boltzmann equation using the relativistic scattering formalism used in [@Rephaeli1995] and [@Phillips1995].
By using a combination of adaptive Gauss-Legendre quadrature and a modified Romberg integration routine, the expressions for the scattered radiation field due to the thermal SZ effect can be solved relatively quickly. 
For example, evaluating the tSZ signal at 1000 frequency points takes on the order of 10 ms.
The kSZ effect requires solving a single integral, which is done solely using an adaptive Gauss-Legendre quadrature.
Because the tSZ and kSZ integrals are solved separately, the cross-terms are added by usng the expansion given by [@Nozawa1998].

In addition to single-pointing signals, `MockSZ` can generate simple spatially-extended models of galaxy clusters, such as the isothermal-$\beta$ model.

Also, `MockSZ` can simulate non-thermal SZ signals from electron populations described by relativistic powerlaws.

# Statement of need
For designing and characterisation of astronomical sub-mm spectrometers, it is important to simulate observational forecasts of signals of interest, taking into account various effects such as atmospheric fluctuations, detector noise, and telescope efficiencies.
`MockSZ` is designed to serve as a source backend for forecasting software such as `TiEMPO` [@Huijten2022]. 

`SZpack` [@Chluba2012], a commonly used SZ simulator, can calculate the SZ signal in multiple different ways, the fastest being expansions of the collisionless Boltzmann equation in appropriate basis functions, such as described in [@Chluba2013] and [@Nozawa2006].
Additionally, `SZpack` can calculate the SZ distortion by directly integrating the Boltzmann equation, but these approaches are rather slow.
`MockSZ` provides an interesting alternative to `SZpack`, working on a different integral equation. The equations solved by `MockSZ` have the advantage of being easy to interpret physically.

# Availability
`MockSZ` is available on [GitHub](https://github.com/arendMoerman/MockSZ) and is available for Linux, MacOS, and Windows.
Software documentation and instructions regarding installation, contributing and issue tracking can be found in the [documentation](https://arendMoerman.github.io/MockSZ/).

# Acknowledgements
This work is supported by the European Union (ERC Consolidator Grant No. 101043486 TIFUUN). 
Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. 
Neither the European Union nor the granting authority can be held responsible for them.

# References
