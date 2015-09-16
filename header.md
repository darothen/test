---
title: Metamodeling of Droplet Activation for Global Climate Models
authors:
    - name: Daniel Rothenberg
      affiliation: 1
      email: darothen@mit.edu
    - name: Chien Wang
      affiliation: 1
      email: wangc@mit.edu
affiliations:
    - id: 1
      address: Massachusetts Institute of Technology, Department of Earth, Planetary, and Atmospheric Sciences
date: May 27, 2015
abstract: 
    The nucleation of clouds droplets from the ambient aerosol is a critical physical process which must be resolved for global models to faithfully resolve aerosol-cloud interactions and aerosol indirect effects on climate. In order to better resolve droplet nucleation from a complex, multi-modal and multi-component aerosol population within the context of a global model, a new metamodeling framework is applied to derive an efficient and accurate activation parameterization. The framework applies polynomial chaos expansion to a detailed parcel model in order to derive an emulator which maps thermodynamic and aerosol parameters to the supersaturation maximum achieved in an adiabatically ascending parcel and can be used to diagnose droplet number from a single lognormal aerosol mode. The emulator requires much less computational time to build, store, and evaluate than a high-dimensional lookup table. Compared to large sample sets from the detailed parcel model, the relative error in the predicted supersaturation maximum and droplet number computed with the best emulator are $-0.6\% \pm 9.9\%$<!-- LARS order 5--> and $0.8\% \pm 17.8\%$<!-- LARS order 5--> (one standard deviation). On average, the emulators constructed here are as accurate and between 10 and 17 times faster than a leading physically-based activation parameterization. Because the underlying parcel model being emulated resolves size-dependent droplet growth factors, the emulator captures kinetic limitations on activation. The results discussed in this work suggest that this metamodeling framework can be extended to accurately account for the detailed activation of a complex aerosol population in an arbitrary coupled global aerosol-climate model.
numbersections: True
graphicspath: figures/
biblio-files: /Users/daniel/Desktop/mendeley_bib.bib
biblio-style: ametsoc2014.bst
header-includes:
    - \usepackage{algorithm}
    - \usepackage{algpseudocode}
    - \usepackage{multirow}
...

<!--
DATA / FIGURE SOURCES:
    1. /Users/daniel/Dropbox (MIT)/Research/slides/sep 5 2014 - single mode (SM, NS2003) analyses; tables of error statistics, LHS big overview plots)
    2. /Users/daniel/workspace/Research/scripts/pce_comparison/ghan_plots
-->
