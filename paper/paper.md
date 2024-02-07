---
title: "Generating mock Sunyaev-Zel'dovich maps with MockSZ"
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
The Sunyaev-Zel'dovich (SZ) effect is an interaction between cosmic microwave background (CMB) photons and the hot, ionised gas residing in clusters of galaxies [@Sunyaev1970].

`MockSZ` is a Python package for simulating on-sky SZ maps of multiple flavours. 
It numerically integrates the scattering and redistribution integrals given in [@Rephaeli1995] and [@Birkinshaw1999], therefore taking into account the relativistic contributions to the total effect. 
In addition to the tSZ and kSZ effects discussed earlier, `MockSZ` can also simulate non-thermal electron distributions by using a relativistic powerlaw

The 

# Statement of need
For designing 

# Availability
`MockSZ` is available on [GitHub](https://github.com/PyPO-dev/PyPO) and is available for Linux, MacOS, and Windows.
Software documentation and instructions regarding installation, contributing and issue tracking can be found in the [documentation](https://pypo-dev.github.io/PyPO/).

# Acknowledgements
This work is supported by the European Union (ERC Consolidator Grant No. 101043486 TIFUUN). 
Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. 
Neither the European Union nor the granting authority can be held responsible for them.

# References
