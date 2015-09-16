# Methodology {#sec:methodology}

## Parcel Model {#sec:methodology-parcel-model}

Adiabatic parcel models are commonly used to study droplet activation and its sensitivity to factors such as environmental conditions and ambient aerosol properties. For this work, a novel parcel model based on previous studies [@Leaitch1986; @Nenes2001; @Seinfeld2006] was designed and implemented to accommodate diverse, chemically heterogeneous, polydisperse aerosol populations. The model simulates droplet growth on the initial aerosol population due to condensation within a constant-speed adiabatic updraft.

Although an arbitrary aerosol size distribution function can be supplied as an input to the model, for the purposes of this study the initial aerosol distribution is assumed to be lognormal and described by the equation

\begin{equation}
n_N(r) = \frac{d N}{d \ln r} = \frac{N_t}{\sqrt{2\pi}\ln \sigma_g}\exp\left(-\frac{\ln^2(r/\mu_g)}{2 \ln^2 \sigma_g}\right) 
\end{equation}

where the parameter set \((N_t, \mu_g, \sigma_g)\) correspond, respectively, to the total aerosol number concentration, the geometric mean radius, and the geometric standard deviation of the distribution. Within the model, this distribution is discretized into 200 size bins equally spaced over the logarithm of particle radius ($r$) and covers the size range \((\text{min}\left[0.1\,\si{nm},\, \mu_g/10\sigma_g\right], \mu_g\times10\sigma_g)\). The mean radius in each bin grows due to condensation so that the activation of wetted aerosol into droplets is calculated in a Lagrangian sense. To relate size-dependent droplet growth to its embedded aerosol's chemical composition, each bin is prescribed a hygroscopicity following $\kappa$-Köhler theory [@Petters2007]. 

To simulate droplet activation, the parcel model first computes an equilibrium wet-size distribution from the given initial aerosol population and initial environmental temperature, pressure, and relative humidity. Then, a set of conservation equations which describe the evolution of the parcel temperature, supersaturation, liquid/vapor water content, and pressure are integrated forward in time using a solver suitable for stiff systems [VODE; @Brown1989]. The complete system of equations and further details on the parcel model can be found in Appendix \ref{app:parcel-model-description}.

## Polynomial Chaos Expansion {#sec:methodology-pce}

We construct an emulator of the parcel model in order to asses droplet activation by applying the probabilistic collocation method [PCM; @Tatang1997]. The PCM maps an output from the parcel model to a set of input parameters by building a response surface using a polynomial chaos expansion. The polynomial which results from this process is optimally a computationally-efficient, high-fidelity reproduction of the detailed parcel model simulation. Although often used for conducting global sensitivity analyses [@Pan1997; @Calbo1998; @Mayer2000; @Lucas2005; @Anttila2007] chaos expansion-based emulators have also been used to build deterministic parameterizations [@Cohen2011]. To apply and build the chaos expansions discussed here, the open-source Design Analysis Kit for Optimization and Terascale Applications [DAKOTA; @Dakota2014] version 6.1 was used, which automates the sampling of the PCM collocation points and the computation of the coefficients of the polynomial chaos expansion given a user-generated interface to a numerical model (the parcel model described in Section \ref{sec:methodology-parcel-model}) and a description of the inputs/outputs to and from that interface.

A review of the theoretical basis of polynomial chaos expansion and its potential applications is provided by @Sudret2008; here, we highlight the important details of the technique as applied via PCM for the benefit of the reader. PCM is a non-intrusive polynomial chaos expansion technique; rather than require complex, significant modifications to the model being emulated,  PCM instead considers the model to be a black box and constructs a map from an input parameter space to the model output parameter space. To accomplish this, PCM re-casts the input parameters to a model as a set of $M$ independent random variables, $\mathbf{X} = {X_1,\dots,X_M}$, each with an associated probability density function. For each input in $\mathbf{X}$, the associated PDF is used as a weighting function to derive an orthogonal polynomial which adds to the bases for the polynomial chaos expansion, $\phi_j$. Using a finite number of these bases, the chaos expansion for a given model response, $R$, is then

\begin{equation}
R \approx \sum\limits_{j=0}^{P} \alpha_j\phi_j(\mathbf{X}) \label{eq:pce_def}
\end{equation}

The complete basis of polynomials up to a fixed total-order $p$ is retained in the expansions computed here. For such a total-order expansion, equation \eqref{eq:pce_def} has $N_t = P + 1 = (M + p)!/(M!p!)$ terms as it contains each of the $p+1$ orthogonal basis polynomials for each input parameter. PCM provides an experimental design for determining the coefficients of the expansion, $\alpha_j$, by evaluating the model response for a set of $N_s$ total input parameter sets, ${\mathbf{X}^1,\dots,\mathbf{X}^{N_s}}$, corresponding to the roots of $\phi_j$ and solving a regression problem

\begin{equation}
\boldsymbol{\Phi}\boldsymbol{\alpha} = \mathbf{R} \label{eq:pce_regress}
\end{equation}

where $\mathbf{R}$ is the vector of model responses, $\boldsymbol{\alpha}$ is the vector of expansion coefficiens, and the matrix $\boldsymbol{\Phi}$ contains rows for each of the polynomial terms $\phi_j$ evaluated for a given input parameter set $\mathbf{X}^j$.

A practical consideration in applying the PCM to a particular problem is what subset of the $N_s$ potential points to use in solving for the coefficients. In general, there exists a full factorial design of size $N_s = (p+1)^M$ available for use (all the roots of the orthogonal basis polynomials for all inputs). However, for even moderately-sized $p$ and $M$, the number of potential model evaluations grows very rapidly. In our application we choose a subset of $N_s'$ parameter sets when applying the PCM by using two rules of thumb:  

1. Choose parameter sets with roots closest to the origin [@Sudret2008];
2. Cross-validate the regression result using $3 N_t$ parameter sets chosen according to rule (1). 

These rules will always produce an over-determined system for Equation \eqref{eq:pce_regress}. The accuracy of the resulting emulators derived in this study were not sensitive to increasing $N_s'$, and $3 N_t$ does not produce an excessive number of required parcel model simulations. Three different techniques were tested for solving this system: typical linear regression by ordinary least squares (OLS), Least Angle Regression [LARS; @Efron2004], and Least Absolute Shrinkage and Selection Operator [LASSO; @Tibshirani2011]. Both LARS and LASSO involve computing 

\begin{equation}
\boldsymbol{\alpha} = \text{arg min}\,||\boldsymbol{\Phi}\boldsymbol{\alpha} - \mathbf{R}||^2_{l_2}\quad \text{such that}\,||\boldsymbol{\alpha}||_{l_1} \leq \tau
\end{equation}

in an iterative, greedy fashion with the potential to yield sparse solutions with some coefficients $\alpha_j = 0$. This would be desirable for high-order chaos expansions for many input parameters, as it would reduce the number of coefficients necessary to save for re-using the expansion as an emulator. 

## Emulation of Parcel Model {#sec:methodology-emulation}

The PCM was applied to emulate the activation of a single, lognormal aerosol mode embedded in a constant-speed adiabatic updraft as simulated by the parcel model described in Section \ref{sec:methodology-parcel-model}. Specifically, the model was used to predict the base-10 logarithm of the maximum supersaturation $S_\text{max}$ given aerosol of different lognormal size distributions and hygroscopicities, for different environmental and thermodynamic conditions. The mechanics of the PCM - and polynomial chaos expansion more generally - permit the use of arbitrary PDFs to describe the input parameters over their physically relevant values.  In this application, uniform distributions were chosen to emphasize that the derived chaos expansion should peform equally well anywhere within the input parameter space. Each uniform distribution is defined by minimum and maximum permissible bounds for each input parameter, $a$ and $b$, such that its probability distribution is just given as $f(x) = 1/(b-a)$ for $a < x < b$. The complete list of input parameters and their bounds is summarized in Table \ref{table:single_mode_params}. 

In order to utilize the PCM, the uniform distribution for each parameter must be rescaled to the range $[-1, 1]$. This produces a new set of random variables for each parameter $X_i$:

\begin{equation}
Z_i = \frac{2(X_i - a_i)}{b_i - a_i} - 1 \label{eq:uni_to_uni}
\end{equation}

The orthogonal polynomials used in the basis of the chaos expansion which correspond to a uniform PDF over the interval $[-1, 1]$ are the canonical Legendre polynomials which follow the three-term recurrence relation:

\begin{align}
    P_0(Z) &= 1 \\
    P_1(Z) &= Z \\
    P_{n+1}(Z) &= \frac{(2n+1)ZP_n(Z) - nP_{n-1}(Z)}{n+1}
\end{align}

The roots of these Legendre polynomials can be inverted using Equation \eqref{eq:uni_to_uni} to determine values in the original, physical parameter space to use in sampling the parcel model. 

The bounds for the physical parameters supplied to the PCM were chosen in order to characterize activation near cloud base (Table \ref{table:single_mode_params}). The logarithm of several variables (aerosol number concentration, aerosol geometric mean radius, and updraft velocity) is used because the supersaturation maximum computed by the parcel model is sensitive to changes in these parameters over several orders of magnitude. Updraft velocity is permitted to vary between $0.01$ and $10.0$ \si{\metre\per\second}; over this range (which covers a spectrum from weakly convecting, stratiform clouds to strong, deeply convecting ones) and the range of aerosol number concentration (which includes clean and very polluted regimes), activated fraction can range from virtually nothing to complete activation of the entire aerosol population. Furthermore, the aerosol mode geometric mean radius $\mu_g$ spans a variety of smaller Aitken-type modes to large, coarse aerosol modes and potentially giant CCN. 

Many parameterizations of droplet nucleation diagnose activation directly by applying equilibrium Köhler theory. To do this, the maximum supersaturation achieved by a cloudy parcel is used as a threshold; particles with a Köhler-predicted critical supersaturation lower than this maximum environmental supersaturation are considered to be activated. However, physically, for a droplet to activate it must grow beyond a critical size corresponding to this critical supersaturation. Due to kinetic limitations on droplet growth, this may not be realized for droplets growing on very large CCN [@Nenes2001]. @Ghan2011 suggests particles with radius larger than 0.1 \si{\micro\metre} or those whose critical supersaturations are close to the environmental maximum supersaturation are likely to suffer from this effect. By directly considering a detailed parcel model, the emulators constructed here consider kinetic limitations on droplet growth and their feedback on the evolving parcel supersaturation. In existing, physically-based parameterizations in the literature, an instantaneous growth-rate assumption must be applied. This assumption causes parameterizations to under-predict supersaturation maximum because instaneous growth will tend to condense water from the vapor phase too quickly and release surplus latent heat which supresses the increase of the supersaturation. 

Because the computed supersaturation maximum in a parcel model activation simulation can also vary over several orders of magnitude, we use $\log_{10} S_\text{max}$ as the response function emulated by the PCM. However, in order to apply this transform to the response function, it must be assumed that the cloudy parcel always supersaturates with respect to water vapor, i.e. $S_\text{max} > 0$. To ensure this, all simulations performed during sampling by the PCM start with an aerosol population equilibrated to 100\% relative humidity and an initial environmental supersaturation of $0$. In very clean conditions with small aerosol particles and weak updraft speeds, this could tend to over-predict the maximum supersaturation (which may not even reach saturation for a sufficiently dry parcel sampled from the boundary layer). Many existing parameterizations in the literature implicitly make this same assumption by representing the aerosol population with respect to a coordinate derived from the critical supersaturation for a given size; in this case the integral over the size distribution spans $0 \leq S \leq S_\text{max}$, and thus considers the same situation with respect to the growth of the nascent droplet population. 

For the 8-parameter input space governing single-mode activation considered here, the 3rd- and 4th-order chaos expansions produced by the PCM have 165 and 495 terms, respectively. The number of terms is equivalent to the number of coefficients one must store in order to re-use a given chaos expansion. This small memory footprint afford chaos expansions a huge advantage over similar parameterizations based on detailed look-up tables. An isotropic look-up table with $M$ parameters and $n$ sample points for each parameter would require $n^M$ values to be stored - a value which for even small numbers of parameters can be several orders of magnitude larger than even a high-order chaos expansion. A more detailed description of how the chaos expansions are saved and later evaluated is given in Appendix \ref{app:pce-emulator}.

From the value of $S_\text{max}$ produced by the parcel model emulator for a given input parameter set, the number concentration of cloud droplest activated, $N_\text{act}$, can be obtained by integrating over the original lognormal aerosol size distribution re-expressed as a function of critical supersaturation rather than droplet radius [@Ghan2011], yielding the expression

\begin{equation}
N_\text{act} = \frac{N}{2}\left(1 - \text{erf}\left[ 2 \ln \left(\frac{S_m}{S_\text{max}}\right)\big/(3\sqrt{2}\ln\sigma_g) \right]\right)
\end{equation}

where $S_m$ is the critical superaturation for the geometric mean radius, $\mu_g$. 

