
# Metamodeling of Droplet Activation for Global Climate Models

**Authors**: Daniel Rothenbeg and Chien Wang

In the archived repository, only the `zip` files for the experimental data are tracked. They expand to rather large archives (~600 mb across all three!) which aren't remotely all necessary for studying; all the sampling experiments are summarized in the `SM_sampling_[results|tidy].csv` file at the top of each data folder. Furthermore, the attributes necessary to evaluate each chaos expansion (the coefficient vector and order matrix) are archived in the top-level HDF5 file, `pcm_param.h5`, which is read by the provided Python code (`pcm_param.py`).

## Abstract

> The nucleation of clouds droplets from the ambient aerosol is a key critical physical process which must be resolved for global models faithfully resolve aerosol-cloud interactions and aerosol indirect effects on climate. In order to better resolve droplet nucleation from a complex, multi-modal and multi-component aerosol population within the context of a global model, a new metamodeling framework is applied to derive an efficient and accurate activation parameterization. The framework applies polynomial chaos expansion to a detailed parcel model in order to derive an emulator which maps thermodynamic and aerosol parameters to the supersaturation maximum achieved in an adiabatically ascending parcel and can be used to diagnose droplet number from a single lognormal aerosol mode. The emulator requires much less computational time to build, store, and evaluate than a high-dimensional lookup table. Compared to large sample sets from the detailed parcel model, the relative error in the predicted supersaturation maximum and droplet number computed with the best emulator are $-0.6\% \pm 9.9\%$<!-- LARS order 5--> and $0.8\% \pm 17.8\%$<!-- LARS order 5--> (one standard deviation). On average, the emulators constructed here are as accurate and between 10 and 17 times faster than a leading physically-based activation parameterization. Because the underlying parcel model being emulated resolves size-dependent droplet growth factors, the emulator captures kinetic limitations on activation. The results discussed in this work suggest that this metamodeling framework can be extended to accurately account for the detailed activation of a complex aerosol population in an arbitrary coupled global aerosol-climate model.
> 
> ## Notes
> 
> 1. Some hand-tweaking is necessary for the JAS submission after generation here:
>   - [ ] rename relative figure paths to just the leaf (figure name)
>   - [ ] something else?