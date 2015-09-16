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
# Introduction 

Interactions between aerosol and clouds yield one of the largest sources of uncertainty in understanding climate and future climate change on regional and global scales [@Boucher2013]. Within the Earth's atmosphere, homogeneous liquid water droplet formation is not thermodynamically favorable [@Pruppacher1997]; instead, the pathway to nucleating cloud droplets is aided by the presence of ambient aerosol, a subset of which possess physical and chemical characteristics which allow them to serve as cloud condensation nuclei (CCN). These CCN provide a linkage between the physiochemical processes of atmospheric particles and cloud microphysics.

Changes in the background aerosol population can directly affect the properties of a nascent cloud droplet population. For instance, holding liquid water content constant, an increase in the number of CCN would tend to increase the total cloud droplet number concentration (the "Twomey" effect)  while necessarily reducing the average size of the droplets [@Twomey1974]. Such a change could enhance a cloud's albedo, an effect which could be further amplified through microphysical feedbacks since smaller droplets impede the production of drizzle and thus lengthen cloud lifetime [@Albrecht1989]. Mechanisms whereby aerosol influence the properties of clouds (and ultimately climate) are generally known as "aerosol indirect effects" [@Haywood2000; @Lohmann2005] and provide a path for changes in the ambient aerosol to produce cascading effects up to progressively larger scales of atmospheric motion [e.g. @Wang2005; @Ekman2011; @Morrison2011; @Tao2011; @Fan2012; @Altaratz2014].

Aerosol indirect effects can either warm or cool the climate, but they all fundamentally depend on a subset of the ambient aerosols which function as CCN. The theory describing the dependency of cloud droplet nucleation (also known as aerosol activation) on CCN availability and ambient aerosol has been rigorously developed using adiabatic and entraining parcel theory [@Seinfeld2006; @Pruppacher1997], and depends on details of the heterogeneous chemical composition, number, size distribution(s) and mixing state of the background aerosol [@Mcfiggans2006] as well as local meteorology [@Morales2010]. Under polluted conditions, effects relating to chemical composition could produce a climatic effect as large as the basic "Twomey" effect [@Nenes2002a; @Lance2004]. 

The development of activation parameterizations was pioneered by @Twomey1959 and @Squires1961, who derived a relationship between the number of activated particles and the environmental supersaturation based on an aerosol size distribution approximated by a power law. @Ghan2011 presented a thorough overview of subsequent developments over the past five decades and an inter-comparison of several modern parameterizations. However, there is still an active effort to improve these parameterizations, as they are increasingly called upon to mediate between ever more complex aerosol models and the climate models to which they are coupled. For instance, the parameterization initially developed by @Nenes2003 has seen continuous development, including modifications to handle condensation onto insoluble but wettable particles using adsorption-activation theory [@Kumar2009], environmental entrainment [@Barahona2007], and numerical improvement of the population-splitting technique [@Barahona2010; @MoralesBetancourt2014b]. Similarly, @Ghan2011 extended the parameterization of @Abdul-Razzak2000 to account for non-unity values of the accommodation coefficient $a_c$. Beyond idealized testing and droplet closure studies [@Meskhidze2005; @Fountoukis2007], these modern parameterizations have been implemented in coupled climate-aerosol models such as the Community Earth System Model (CESM) to predict online cloud droplet number concentrations, where they have been shown to correct biases in global average cloud droplet number concentrations and improve agreement with cloud properties measured from satellite-born instruments [@Gantt2014]. Furthermore, adjoints of these parameterizations have been derived and coupled to chemical transport and global models in order to study the sensitivities of cloud droplet number to aerosol, chemical and microphysical factors [@Karydis2012; @Moore2013].

Additionally, following the original integral/geometric approach by @Twomey1959, analytical representations of supersaturation evolution from adiabatic parcel theory have been progressively generalized to relate aerosol distributions to activation kinetics [@Cohard1998; @Khvorostyanov2006; @Khvorostyanov2008; @Shipway2010; @Shipway2015]. Although fundamentally analytical parameterizations, schemes of this class typically rely on expensive numerical operations, such as in the evaluation of hypergeometric functions and iterative loops.

While most of these recent efforts towards improving activation parameterizations have focused on building highly generalized, "physically-based" tools, there is still an application for other parameterization approaches. @Saleeby2004 parameterized droplet nucleation for a cloud-resolving model, the Regional Atmospheric Modeling System (RAMS), by constructing a four-dimensional lookup table based on temperature, vertical velocity, aerosol number concentration, and the median radius of a lognormal aerosol mode with assumed chemical composition. A fifth dimension representing chemical composition via aerosol hygroscopicity (following $\kappa$-Köhler theory [@Petters2007]) was added to the lookup table [@Ward2010] and later generalized to aerosol soluble fraction [@Saleeby2013]. Constructing lookup tables of detailed parcel model results can be considered a form of model emulation combining a cache of known results and local polynomial (linear) approximation. 

As the degrees of freedom and number of parameters describing a given aerosol population in a model increase, the burden of saving enough known points to interpolate through the parameter space via lookup table to some reasonable accuracy increases algebraically. For instance, the CESM features a modal aerosol population with three predefined, internally mixed lognormal modes, each with a fixed geometric standard deviation [@Liu2012]. Each mode is uniquely described by two moments (total number and total mass concentration) and the chemical composition of the mode by a single prognostic hygroscopicity term. Thus, the entire aerosol population has $N=9$ degrees of freedom - too many to build a lookup table of activation statistics. Physically-based parameterizations were designed to accommodate these sorts of arbitrary mixtures of aerosol, but have a tendency to systematically underestimate activated fractions and subsequently cloud droplet number [@Simpson2014]. This is due to the parameterizations' use of a set of assumptions which become increasingly likely to be violated as the aerosol population becomes more complex, specifically (1) that proto-droplets grow in equilibrium with environmental changes in relative humidity, and (2) that there are no kinetic or inertial limitations to droplet growth. The presence of giant CCN [@Barahona2010] and weak updrafts or excessively polluted conditions [@Nenes2001] exacerbates this problem. 

The goal of this study is to apply a surrogate modeling or emulation technique commonly used in the uncertainty quantification literature to a detailed parcel model capable of describing aerosol activation; this yields an efficient parameterization optimized for the high-dimensional parameter space effecting droplet nucleation in coupled aerosol-climate model. In essence, employing the derived emulator as an activation parameterization would be akin to directly coupling a detailed parcel model to a global model. Such a parameterization would be directly physically-based, but rely on fewer assumptions which affect the condensational growth of aerosol into CCN. However, it would also incorporate the efficiency of a lookup table, since the emulator would be designed to require a scarce amount of cached information and to be computationally cheap to evaluate. Additionally, it would improve upon the framework of lookup table and be extensible to a very high-dimensional parameter space and thus be compatible with aerosol-climate models of increasing complexity. 

Section \ref{sec:methodology} describes the parcel model and the probabilistic collocation method (PCM) used to build its emulator. Section \ref{sec:results} presents results from applying the PCM to build a parcel model emulator designed to simulate the activation of a single lognormal aerosol mode under a wide variety of background environments and compares the new emulator to existing activation parameterizations. Section \ref{sec:conclusion} motivates an extension of the technique to an emulator suitable of mediating aerosol activation in a coupled aerosol-climate model. 

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

# Results {#sec:results}

## Evaluation of Emulators {#sec:results-evaluation}

In order to assess the performance of the emulator, two sets of $n=10,000$ samples were drawn using maximin Latin Hypercube Sampling from the parameter space defined in Table \ref{table:single_mode_params}. This randomized design ensured that representative, equal numbers of samples were drawn from across the multi-dimensional parameter space. In the first set, variables whose logarithms were used to build the emulator were sampled in logarithmic space; in the second set, these variables were transformed back to their original values (e.g. from $\log_{10} N$ to $N$) before the sample was constructed. The two independent sets were blended together to assess the emulator. This helps ensure that both very high and very low values of the log-transformed are thoroughly represented within the sample. The set of sample parameter sets were run through all the derived chaos expansions of all orders, as well as the detailed parcel model as a reference benchmark for activation dynamics.

Figure \ref{fig:single_mode_LHS} illustrates the performance of a 4th-order expansion whose coefficients were derived using ordinary least squares. The large range of initial temperatures, pressures, aerosol populations, and updraft speeds sampled here leads to a very large range of supersaturation maxima achieved by the ascending parcel. Weaker updraft speeds are generally associated with lower maximum supersaturations and corresponding to lower aerosol activated fraction; the opposite is true when strong updrafts are present, although there are some cases where a strong updraft activates a small fraction of aerosol. This typically occurs when initial aerosol size distribution is shifted towards larger radii and under polluted conditions with aerosol number concentrations greater than $3,000$ \si{\per\cubic\centi\metre}. However, over the large parameter space sampled, the chaos expansion accurately reproduces the parcel model's determination of $S_\text{max}$ and corresponding activated fraction. 

For parameter sets leading to a large activated fraction of 0.8-1.0, the relative error of the chaos expansion (compared to the parcel model) rarely exceeds 5\%, and on average (for all activated fractions) is 5.7\%. While the mean relative error in each activated fraction decile is close to 0, the standard deviation in the relative error tends to increase for the lower ones; the standard deviation in relative error decreases from 19.7% for activated fractions ranging from 0.1-0.2 to 3.5\% for those ranging from 0.8-0.9. This suggests that there is a non-linear component in the mapping from the input parameter space to the emulated maximum supersaturation and diagnosed droplet number concentration which is prevalent in the weak droplet activation regime; the predicted activated fraction is more sensitive to small changes in the input parameters in this regime than in others.

Increasing the order of the chaos expansion tends to improve the accuracy of the predicted $S_\text{max}$, as recorded in Table \ref{table:smax_stats}. However, there is not much difference between the methods used to compute the coefficients of the expansion beyond expansion order. For example, for the 4th- and 5th-order expansions, the expansions perform equally well regardless of what method (OLS, LARS, or LASSO) was used to compute the coefficients when considering the mean and spread of the relative error to the parcel model reference simulations. In all cases, the chaos expansions produce very large r$^2$ values and normalized fractional root mean square errors ($\text{RMSE}/\sum\limits_{i=1}^n(X_i^2/n)$) which decrease as the order of the expansion increases.

These same statistics, computed for the diagnosed droplet number concentration given the predicted supersaturation maximum, are summarized in Table \ref{table:Nd_stats}. Here, the trend is similar to before; increasing the order of the expansion tends to improve the accuracy of the diagnosed number concentration in terms of mean relative error, and also tends to decrease spread around that value. Third-order expansions tend to produce more accurate results with respect to the mean relative error, but this is overshadowed by the fact that there is far more variance in their predicted values as indicated by the standard deviation of their relative errors, which are almost twice as large as those of the higher-order expansions. 

## Comparison with Other Parameterizations {#sec:results-comparison}

We compare the performance of the chaos expansion-based emulators to two existing parameterizations from the literature. The scheme by @Abdul-Razzak2000 (ARG) - which is widely used in global models - utilizes a psuedo-analytical solution to an integro-differential equation derived from the adiabatic parcel system with embedded aerosol growing via condensation. This is in contrast to the scheme by @MoralesBetancourt2014b (MBN), which instead utilizes an iterative scheme to separate the aerosol population into subsets whose growth is inertially limited or not, and uses this information to derive a maximum supersaturation for a given parcel system. The MBN scheme is generally more expensive to evaluate than the ARG scheme due to its iterative nature, but is often more accurate owing to its consideration of the potentially important effect of large albeit unactivated aerosol particles [@Simpson2014]. In contrast with these schemes, our emulators simulate the activation process based on the explicit numerical solution obtained from a detailed parcel model which is similar to the one used to build and evaluate the MBN scheme [e.g. @Nenes2003] 

The parameter sets used in Section \ref{sec:results-evaluation} were also used to compute droplet activation with the ARG and MBN parameterizations. The relative errors between the supersaturation maximum and droplet number predicted by these schemes and the chaos expansions compared to the parcel model are illustrated in Figure \ref{fig:mre_boxplots}. Both the ARG and MBN parameterizations are more accurate than the 2nd and 3rd order chaos expansions. The ARG scheme tends to underpredict the maximum supersaturation, which is consistent with previous investigations into its performance [@Abdul-Razzak2000; @Ghan2011; @Simpson2014]. This tends to produce a bias towards under-prediction of droplet number. The MBN scheme tends to yield more accurate predictions of both maximum supersaturation and droplet number. However, both schemes are out-performed by the 4th and 5th order chaos expansions, both on average and in terms of the variance of the predictions; for instance, the OLS-derived 4th and 5th order expansions yield relative error in predicted $S_\text{max}$ with a mean and standard deviation of $4.2\% \pm 14.4\%$ and $-0.32\% \pm 10.4\%$, whereas the ARG and MBN schemes yield $-13.7\% \pm 18.2\%$ and $-2.4\% \pm 15.8\%$, respectively. This is larger than other studies have reported, but we explore a much larger parameter space in our sampling for the purposes of deriving the chaos expansion.

As a consequence of tending to slightly under-predict $S_\text{max}$, both the ARG and MBN schemes underpredict the number of activated droplets in the framework considered here. The mean relative error in droplet number predicted by the ARG and MBN schemes for the samples here are $-9.6\% \pm 23.4\%$ and $-4.9\% \pm 16.8\%$, respectively. All of the chaos expansions outperform the mean relative error of the ARG scheme, and those of order $p \geq 3$ do so with less variance. 

In addition to producing low mean relative error in predicted $S_\text{max}$ and droplet number activated, the chaos expansions also reproduce the dependence of activation dynamics on aerosol physical properties and updraft speed, as illustrated in Figure \ref{fig:single_var_sensitivity}. In response to increasing aerosol number concentration, the $S_\text{max}$ reached by an ascending parcel tends to decrease as enhanced competition by aerosol for water vapor produces a larger source of latent heat release and warming, which limits the production of supersaturation; overall this increases the droplet number concentration as the aerosol activated fraction only decreases by a factor of 4 when the total number of initial aerosol increases by an order of magnitude (Figure \ref{fig:single_var_sensitivity}a-b). Shifting the aerosol population to larger sizes (Figure \ref{fig:single_var_sensitivity}c-d) produces a similar effect in inhibiting the increase in a parcel's supersaturation; however, Köhler theory predicts that these larger particles will more easily activate, which offsets the increase in $S_\text{max}$ and yields larger droplet number concentrations. A similar effect occurs as aerosol hygroscopicity increases (Figure \ref{fig:single_var_sensitivity}e-f). 

The chaos expansions, as well as both the ARG and MBN schemes, capture these subtleties of activation dynamics as well as the detailed parcel model. More importantly, the expansions reproduce the sensitivity of activation to updraft speed (Figure \ref{fig:single_var_sensitivity}g-h) which is the most important factor controlling $S_\text{max}$ and setting the droplet number. At the largest updraft speeds of a few meters per second - indicative of deep, vigorous convection - the MBN scheme outperforms both the chaos expansions and the ARG scheme. However, for the aerosol population considered in Figure \ref{fig:single_var_sensitivity}g-h (with $N = 1000\,\si{\per\cubic\cm}$, $\mu = 0.05\,\si{\micro\metre}$, and $\sigma=2.0$), the relative error in predicted $S_\text{max}$ by the chaos expansions at high updraft speeds does not substantially affect the diagnosed droplet number concentration, since in this case all but the smallest aerosol particles activate under equilibrium considerations. Note that this $S_\text{max}$-overprediction coupled with an accurate assessment of activated fraction occurs for many different single-mode, lognormal aerosol populations.

Since Figure \ref{fig:single_var_sensitivity} highlights the fact that different schemes potentially perform better in different parts of the parameter space governing droplet activation, we stratified the sampling results based on level of pollution and updraft speed strength and computed activated fraction relative error statistics in each of these bins as shown in Figure \ref{fig:stratified_error}. All of the schemes are accurate in clean and lightly polluted conditions (with aerosol number concentration, $N < 1000\,\si{\per\cubic\cm}$). However, there is a tendency for both the ARG and MBN schemes to underpredict droplet number in heavily polluted conditions ($N > 2500\,\si{\per\cubic\cm}$). 

The 4th-order OLS-derived chaos expansion is plotted in Figure \ref{fig:stratified_error}a as a representative example of the chaos expansions, and it retains its accuracy across the pollution level-updraft strength spectrum. The combination of light updrafts and heavy pollution tends to produce the largest under-prediction in activated droplet number, ranging from 10-30\% for the ARG and MBN schemes. Unsurprisingly, relative error in activated fraction tends to be least sensitive to increasing aerosol number concnetration in the strong updraft regime. In this case, the vigorous updraft produces strong adiabatic cooling which overwhelms latent heat release from condensation as the droplets in the parcel grow, producing a large $S_\text{max}$ and thereby activating a significant fraction of the aerosol.

It should be noted that parts of the parameter space we considered in applying the PCM, evaluating its output, and comparing to existing parameterizations may not be typical of real atmospheric cases. We chose a large parameter space in order to derive the most general emulator possible for this particular single mode aerosol case. In the real world, there should be some correlation between the ambient temperatures, pressures, and updraft speeds used when diagnosing aerosol activation, while we sample these factors as if they were independent from one another. The output from applying the parcel model to non-realistic activation scenarios could tend to inflate the computed relative errors and exacerbate differences between the parcel model and the existing, physically-based parameterizations.  

## Computational Efficiency of Chaos Expansions {#sec:comp-efficiency}

As detailed in Appendix \ref{app:pce-emulator}, evaluating the emulator produced by the chaos expansion requires two sets of straightforward floating-point operations. The first set requires the projection of the the input parameters into the vector space spanned by the basis polynomials of the chaos expansion, which can then be used to evaluate the basis polynomials up the required order. The remaining operations simply multiply these intermediate evaluations together and sum them to evaluate the full expansion. In general, this procedure should lie in between the ARG and MBN schemes in terms of computational complexity. The ARG scheme relies on straightforward floating point operations to derive an estimate for $S_\text{max}$ which at worst involve evaluating a logarithm. However, the MBN scheme requires sets of iterations, each of which necessitates a costly evaluation of the error function and the complementary error function. 

On average, the 2nd and 5th-order OLS-derived chaos expansion was 10-17 times faster than the MBN scheme given the same single mode aerosol population. The exact speedup depended on the background velocity; for weak updraft speeds, the performance of the MBN scheme fared better, although it became much worse for updrafts where $V < 2\,\si{\meter\per\second}$. The ARG scheme was consistently 1-3 times faster than those same chaos expansions. Since the pathway for evaluating either the ARG or chaos expansion schemes do not change depending on the input parameters, their performance was the same regardless of what inputs were provided. 

# Summary and Conclusions {#sec:conclusion}

An efficient parameterization of droplet activation for a single aerosol model under a wide variety of different physico-chemical properties and thermodynamic conditions was developed via statistical emulation of a detailed parcel model using polynomial chaos expansion. The emulators predict the maximum supersaturation achieved by a parcel, which is then used to diagnose activated droplet number using Köhler theory in a similar framework to existing activation parameterizations. The 4th- and 5th-order chaos expansions derived from the detailed parcel model are more accurate on average than two commonly-used, physically-based parameterizations from the literature [@Abdul-Razzak2000; @MoralesBetancourt2014b]. Additionally, the chaos expansions are all at least 10 times faster to evaluate than the MBN scheme and only about twice as expensive as the ARG scheme. A simple algorithm was suggested for evaluating a chaos expansion which requires a minimal amount of data about the expansion (such as the basis polynomials and the coefficients of the expansion terms) to be saved; in this way, the chaos expansions offer a method for extending lookup tables to very high dimensionalities without suffering from exponentially-rising storage costs.

Based on the large set of aerosol properties and thermodynamic conditions we sampled in order to derive and evaluate the chaos expansions, we observed that our emulators particularly outperform the existing schemes in conditions where a light updraft and heavy aerosol pollution (with respect to number concentration) are present. Because the ultimate goal of an activation parameterization is to couple the aerosol physics and chemistry to the cloud microphysics of a global-scale model, this deficiency in the existing parameterizations could be particularly important. Few global models have aerosol-microphysics connections in their deep convection parameterizations, but many source potential cloud droplet formation based on a detailed aerosol activation calculation for their shallow convection and stratiform microphysics schemes. These schemes sometimes artificially restrict the lowest possible updraft speed available for estimating droplet activation, but as a consequence they ensure that weak updrafts make up a large portion of the activation conditions considered during a model run. In regions of the world with heavy anthropogenic aerosol pollution - such as southern and eastern Asia - this provides a recipe for systematically under-predicting droplet number and potentially impacting either a global model's simulated aerosol indirect effect on climate or the modeled aerosol-cloud interaction's sensitivity to changes in anthropogenic aerosol emissions.

Critically, the framework from which the chaos expansions reported here are derived is extendable to the case where a complex, multi-species/multi-modal aerosol population is tracked by a global model; in that case, the number of parameters describing the aerosol size distribution and chemical composition simply increases. Future work will derive chaos expansions emulating activation for a multi-modal aerosol distribution specific to a particular global aerosol-climate model. Additionally, physical processes not considred here can also be introduced into the chaos expansion framework. For instance, entrainment can be incorporated into the parcel model following @Seinfeld2006 and @Barahona2007. Subgrid-scale variability in updraft speeds due to the coarse resolution of global model grids and the distribution of these updrafts can be represented either by a characteristic value [@Morales2010] or by numerical integration over a distribution [@Lohmann1999]. In the latter case, many activation calculations must be performed, incurring a large computational cost. However, the entire integration over a spectrum of droplet speeds could be parameterized in the chaos expansion framework, greatly reducing the cost of this calculation and potentially improving the accuracy of diagnosed cloud droplet number.

As the complexity of global aerosol-climate models increases with respect to the number of aerosol modes and species tracked by the model, there is a pressing need to understand how biases in activation calculations across the high-dimensional parameter spaces defining the aerosol-climate model affect cloud properties and ultimately impact modeled climate. This work highlights a novel way to build efficient, accurate activation schemes for this purpose akin to customized lookup tables which cannot themselves extend to cover the necessary parameters. Employing such schemes should help improve simulated cloud microphysical properties and constrain modeled aerosol indirect effects on climate.


# Acknowledgments

The work in thus study was supported by the National Science Foundation
Graduate Research Fellowship Program under NSF Grant No. 1122374. We thank Steve Ghan (PNNL) and Athansios Nenes (Georgia Tech) for reference implementations of their activation parameterizations.\appendix

# Parcel Model Description {#app:parcel-model-description}

The adiabatic cloud parcel model implemented for this study follows the basic equations of @Pruppacher1997 and adopts the framework used by @Nenes2001 to account for kinetic limitations on droplet growth. Fundamentally, the model integrates a system of coupled ordinary differential equations which describe the thermodynamic evolution of an adiabatically-lifted, non-entraining parcel. In all the simulations described here, we use the Variable-coefficient Ordinary Differential Equation solver [VODE; @Brown1989] to integrate the system forward in time. 

The model tracks the evolution of supersaturation, $S$, with respect to water as

\begin{equation}
\frac{dS}{dt} = \alpha(T, P) V - \gamma(T, P) \frac{dw_c}{dt} \label{eq:ds_dt}
\end{equation}

where $\alpha(T, P) = \frac{gM_wL}{c_pRT^2} - \frac{gM_a}{RT}$ and $\gamma(T, P) = \frac{P M_a}{e_sM_w} + \frac{M_wL^2}{c_pRT^2}$ are functions which are weakly dependent on temperature and pressure [@Leaitch1986], $M_w$ and $M_a$ are the molecular weights of water and air, $L$ is the latent heat of vaporation of water, $c_p$ is the specific heat of dry air at constant pressure, $R$ is the universal gas constant, $g$ is the acceleration due to gravity, $e_s$ is the saturation vapor pressure, and $w_c$ is the liquid cloudwater mass mixing ratio. Equation \eqref{eq:ds_dt} expresses the supersaturation as a balance between production due to adiabatic cooling and loss due to latent heat release. This same framework describes the parcel's change in temperature over time,

\begin{equation}
\frac{dT}{dt} = -\frac{gV}{c_p} - \frac{L}{c_p}\frac{dw_v}{dt} \label{eq:dT_dt}
\end{equation}

where $V$ is the updraft velocity and $w_v$ is water vapor mass mixing ratio. Water mass is conserved as vapor condenses into cloud water,

\begin{equation}
\frac{dw_v}{dt} + \frac{dw_c}{dt} = 0 \label{eq:water_cons}
\end{equation}

Equations \eqref{eq:ds_dt}-\eqref{eq:water_cons} are linked through the growth of the cloud droplet population from the initial aerosol. Given $n$ bin sizes, each associated with a number concentration $N$ and a radius $r$, the change in cloud water can be written as 

\begin{equation}
\frac{dw_c}{dt} = \frac{4\pi\rho_w}{\rho_a}\sum\limits_{i=1}^n N_i r_i^2\frac{dr_i}{dt} \label{eq:dwc_dt}
\end{equation}

where $\rho_w$ and $\rho_a$ denote the density of water and air, respectively. 

The diffusional growth rate for droplets in the $i$th bin is calculated by

\begin{equation}
\frac{dr_i}{dt} = \frac{G}{r_i}\left(S - S_\text{eq}\right) \label{eq:dri_dt}
\end{equation}

where $S$ is the environmental supersaturation, $S_\text{eq}$ is the Köhler-predicted equilibrium supersaturation of the droplet, and $G$ is a growth coefficient which is a function of both the physical and chemical properties of the particle receiving condensate,

\begin{equation}
G = \left(\frac{\rho_wRT}{e_sD'_vM_w} + \frac{L\rho_w[(LM_w/RT) - 1]}{k'_aT}\right)^{-1}
\end{equation}

Non-continuum effects on the diffusivity ($D'_v$) and thermal conductivity ($k'_a$) factors are accounted for with the corrections

\begin{equation}
D'_v = D_v\bigg/\left(1 + \frac{D_v}{a_c r}\sqrt{\frac{2\pi M_w}{RT}}\right) 
\end{equation}

and

\begin{equation}
k'_a = k_a\bigg/\left(1 + \frac{k_a}{a_T r \rho_a c_p}\sqrt{\frac{2\pi M_a}{RT}} \right)
\end{equation}

In these correction terms, the thermal accommodation coefficient, $a_T$, is assumed to be 0.96; the condensation coefficient, $a_C$, is allowed to vary as observations suggest it could take values between $0.1-1.0$ [@Raatikainen2013]. The instantaneous droplet growth rate is further modulated by the difference between the environmental supersaturation, $S$, and the saturation ratio over the surface of the aqueous droplet, $S_\text{eq}$. We treat the droplet-dependent $S_\text{eq}$ following @Petters2007, who employ a single-term $\kappa$ to parameterize particle hygroscopicity; values of $\kappa$ can be derived from laboratory experiments. Under the framework of $\kappa$-Köhler theory the curvature effect term remains the same, while the solute effect term is re-written such that 

\begin{equation}
S_\text{eq} = \frac{r^3 - r_d^3}{r^3 - r_d^3(1-\kappa)}\exp\left(\frac{2M_w \sigma_w}{RT\rho_wr}\right)
\end{equation}

where $r$ and $r_d$ are the droplet radius and the dry radius of its embedded aerosol particle (which is tracked for each initial aerosol size in the model), and $\sigma_w$ is the droplet surface tension, which we take to be independent of the droplet solution composition and described following the recommendation of @Pruppacher1997, $\sigma_w = 0.0761 - 1.55\times 10^{-4}(T - 273)$. A limitation of the this approach for computing $S_\text{eq}$ is that it is not convenient to derive analytical expressions for the critical supersaturation and radius; they must be computed numerically by finding the value $r_\text{crit}$ such that

\begin{equation}
\frac{\partial S_\text{eq}}{\partial r}\bigg|_{r_\text{crit}} = 0
\end{equation}

and then computing $S_\text{crit} = S_\text{eq}(r_\text{crit})$ for a given $\kappa$ and $r_d$. This is accomplished using Brent's method [@Brent1973] and by bounding $r_\text{crit}$ from below with the observation that $r_\text{crit} > r_d$.

Finally, the parcel thermodynamic description is closed by predicting the pressure change within the ascending parcel following the hydrostatic relationship, which can be written using the ideal gas law as 

\begin{equation}
\frac{dP}{dt} = - \frac{gPV}{R_dT_v} \label{eq:dP_dt}
\end{equation}

where $T_v$ is the virtual temperature, which is employed to account for changes in air density due to loss of water vapor to condensate. Equations \eqref{eq:ds_dt}-\eqref{eq:dwc_dt}, \eqref{eq:dP_dt}, and \eqref{eq:dri_dt} applied to each $n$ droplet size bins form a closed system which conserves total water mass.

# Chaos Expansion Emulator Evaluation {#app:pce-emulator}

The PCM as applied here produces two outputs: a $P$-length vector of coefficients $\alpha$ comprised of real values<!--$\alpha \in \mathbb{R}^P$-->, and a $P\times M$ matrix of orthogonal polynomial orders $\Phi$ comprised of integers<!--$\Phi \in \mathbb{Z}^{P\times M}$-->. Each term in matrix $\Phi$ contains a multi-index component for each term in the chaos expansion and indicates the order of the orthogonal polynomial corresponding to term $1 \le j \le M$ for expansion term $0 \leq i \leq P$. For any expansion, $\text{max}[\Phi] = p$, the desired order of the chaos expansion. Algorithm \ref{alg:pce_evaluation} describes the evaluation of a chaos expansion.

\begin{algorithm} \label{alg:pce_evaluation}
    \caption{Psuedo-code for evaluating a polynomial chaos expansion of the form given in Equation \eqref{eq:pce_def}, applied to the computation of $S_\text{max}$}
    \begin{algorithmic}[1] 
        \ForAll {$X_j$}
            \State $Z_j \gets \text{project }X_j$ 
        \EndFor
        \State $\hat{S} \gets 0$
        \For {row $i=0$; $i \le P$}
            \State $\hat{S_i} \gets 1$
            \For {column $j=0$; $j \le M$}
                \State $k \gets \Phi[i, j]$ 
                \State $P_{j,k} \gets \phi_j^\text{k}(Z[j])$
                \State $\hat{S_i} = \hat{S_i}\times P_{j,k}$
            \EndFor 
            \State $\hat{S} = \hat{S} + \alpha[i]\times\hat{S_i}$
        \EndFor
        \State $S_\text{max} \gets 10^{\hat{S}}$
    \end{algorithmic}
\end{algorithm}

Evaluating the chaos expansion involves two parts. First, the input parameters must be projected to conform to the space supported by the PDFs associated with each basis polynomial type, producing a set of parameters $Z_j$. In general, a set of mixed orthogonal polynomials could be used to derive a chaos expansion, but here only Legendre polynomials were used, and each parameter can be projected using Equation \eqref{eq:uni_to_uni}. Second, the polynomial can be evaluated by treating $\Phi$ as a look-up table for the orders of each basis orthogonal polynomial. In practice, the evaluation of these orthogonal polynomials at $Z_j$ for orders up to $k$ can be efficiently precomputed (before the polynomial evaluation loop) by existing orthogonal polynomial libraries [@Gautschi1994].

In general, the only computationally complex part of the chaos expansion evaluation algorithm is the projection from $X_j$ to $Z_j$; given certain basis orthogonal polynomials and their associated PDFs, this procedure could involve numerical integration or otherwise complicated function evaluations.  In the case of simple uniform PDFs, though, the process is acheived entirely by re-scaling the parameters, with little computational overhead. Furthermore, although evaluating a chaos expansion requires looping over each of its terms, each term can be computed independently from one another and efficiently optimized. This is in constrast with an iterative scheme which could involve numerical integration or other costly operations and must be performed in sequence.
        <!--
Summary table of input parameters and ranges of values
-->
\begin{table}[ht]
    \begin{center}
    \caption{Input parameters and bounds used in computing chaos expansion for emulating droplet activation from a single, lognormal aerosol mode embedded in a constant-speed updraft.}
    \label{table:single_mode_params}
    \begin{tabular}[ht]{c | c | c | c | c |}
    Symbol & Name & Units & Bounds \\\hline
    $\log_{10} N$ & Mode Number Concentration & $\log_{10}$ \si{cm\tothe{-3}} & [1, 4] \\ 
    $\log_{10} \mu_g$ & Mode Geometric Mean Radius & $\log_{10}$ \si{\micro\metre} & [-3.0, 1.0] \\
    $\sigma_g$ & Mode Standard Deviation & - & [1.2, 3.0] \\
    $\kappa$ & Mode Hygroscopicity & - & [0.0, 1.2] \\  
    $\log_{10} V$ & Updraft Velocity & $\log_{10}$ \si{\metre\per\second} & [-2.0, 1.0] \\
    $T$ & Air Temperature & K & [240, 310] \\
    $P$ & Air Pressure & Pa & [50000, 105000] \\
    $a_c$ & Accommodation Coefficient & - & [0.1, 1.0]
    \end{tabular}
    \end{center}
\end{table}

<!--
Summary statistics - predicted Smax for all chaos expansions
-->
\begin{table}[ht]    
    \caption{Summary statistics for supersaturation maxima predicted by chaos expansions derived in this study compared to detailed parcel model calculations. The chaos expansions are organized by the method used to derive their coefficients and the expansion order in the two left-most columns. From left-to-right, the reported statistics are the normalized root mean square error (NRMSE), coefficient of determination (r$^2$, mean relative error (MRE), and standard deviation of the mean relative error (MRE std. dev.)}
    \label{table:smax_stats}
    \begin{center}
    \begin{tabular}{cc|cccc}
     & & NRMSE &        r$^2$ &       MRE &    MRE std. dev.\\\hline
    \multirow{4}{*}{LASSO} 
     & 2 & 0.292766 &  0.877600 & -4.963330 &  31.098459 \\
     & 3 & 0.193247 &  0.946670 & -1.012509 &  19.797674 \\
     & 4 & 0.112323 &  0.981983 &  3.941698 &  13.396632 \\
     & 5 & 0.119231 &  0.979699 &  0.772525 &   9.904727 \\\hline
    \multirow{4}{*}{LARS} 
     & 2 & 0.325442 &  0.848752 &  3.055232 &  39.776865 \\
     & 3 & 0.250264 &  0.910559 & -5.433902 &  17.853369 \\
     & 4 & 0.104401 &  0.984435 & -0.517486 &  13.092187 \\
     & 5 & 0.135634 &  0.973729 & -0.572383 &   9.962430 \\\hline
    \multirow{4}{*}{OLS} 
     & 2 & 0.266774 &  0.898368 &  7.956622 &  40.762613 \\
     & 3 & 0.220034 &  0.930861 & -1.232172 &  20.441521 \\
     & 4 & 0.201934 &  0.941768 &  4.243836 &  14.483854 \\
     & 5 & 0.128687 &  0.976351 & -0.259165 &  10.430909 \\\hline
    \end{tabular}

    \end{center}
\end{table}

<!--
Summary statistics - predicted Nderiv vs parcel Neq
-->
\begin{table}[ht]
    \begin{center}    
    \caption{Same as Table \ref{table:smax_stats}, but for predicted droplet number concentration.}
    \label{table:Nd_stats}
    \begin{tabular}{cc|cccc}
     & & NRMSE &        r$^2$ &       MRE &    MRE std. dev.\\\hline
     \multirow{4}{*}{LASSO}
     & 2 & 0.165043 &  0.929439 &  5.697161 &  49.813356 \\
     & 3 & 0.108790 &  0.969342 & -0.930886 &  20.270639 \\
     & 4 & 0.074406 &  0.985659 &  2.545998 &  19.834428 \\
     & 5 & 0.058872 &  0.991022 &  1.441665 &  14.497129 \\\hline
     \multirow{4}{*}{LARS} & 2 & 0.167168 &  0.927611 &  3.524631 &  39.031227 \\
     & 3 & 0.121896 &  0.961510 & -0.030537 &  34.613813 \\
     & 4 & 0.075984 &  0.985044 & -0.280769 &  16.430850 \\
     & 5 & 0.058961 &  0.990995 &  0.884227 &  17.817749 \\\hline
     \multirow{4}{*}{OLS} & 2 & 0.174550 &  0.921076 &  7.640377 &  43.978774 \\
     & 3 & 0.125045 &  0.959496 &  0.380943 &  28.562715 \\
     & 4 & 0.079971 &  0.983434 &  2.556756 &  20.929598 \\
     & 5 & 0.061762 &  0.990119 &  1.295857 &  17.903202 \\\hline
    \end{tabular}

    \end{center}
\end{table}

<!--
Figure 1: Mean relative error boxplots
-->
\begin{figure}[t]
    \noindent
    \includegraphics[width=\textwidth]{figures/1_LHS_SAMPLING/boxplots_multi.pdf} 
    \caption{Boxplots illustrating mean relative error between supersaturation max (a) and droplet number concentration (b) predicted by chaos expansions and parameterizations versus detailed parcel model. The chaos expansions have been grouped by expansion order (x-axis) and method for computing their coefficients (OLS, LARS, and LASSO; hue). ``ARG'' refer to the scheme of \cite{Abdul-Razzak2000}; ``MBN'' refers to the scheme of \cite{MoralesBetancourt2014b}.}
    \label{fig:mre_boxplots}
\end{figure}

<!--
Figure 2: Single mode LHS results, using OLS (p=4) as an example. This
          figure will have two panels:
              a) Smax
              b) Nact   
-->
\begin{figure}[t] 
    \noindent
    \includegraphics[width=\textwidth]{figures/1_LHS_SAMPLING/oneone_OLS_4.pdf} 
    \caption{One-one plot comparing predicted supersaturation maximum (a) and diagnosed equilibrium droplet activated fraction (b) between parcel model and a polynomial chaos expansion of order $p=4$, with coefficients computed using ordinary least squares. Black lines denote factor of $2$ difference between predicted values using parcel model and those computed with the parameterization. Glyph shading denotes updraft velocity, $V$, with corresponding scale on panel (a).} 
    \label{fig:single_mode_LHS}
\end{figure}

<!--
Figure 3: Single-variable sensitivity analysis - mu, N, V, kappa
    /Users/daniel/workspace/Research/scripts/pce_comparison/ghan_plots
-->
\begin{figure}[t]
    \noindent
    \includegraphics[width=\textwidth]{figures/2_SENSITIVITY/1d_sensitivity_multi.pdf}
    \caption{Sensitivity of parameterized and simulated maximum supersaturation (a-d) and activated number fraction (e-h) to changes in mode number concentration, mode geometric mean radius, mode hygroscopicity, and updraft speed with all other parameters held fixed at the values $T = 283\,\si{K}$, $P = 850\,\si{\hecto\pascal}$, $V = 0.5\,\si{\meter\per\second}$, $a_c = 1.0$, $\mu = 0.05\,\si{\micro\meter}$, $\kappa = 0.54$, $N = 1000\,\si{\per\cubic\cm}$, and $\sigma=2.0$.  ``MBN'' and ``ARG'' correspond to the schemes of \cite{MoralesBetancourt2014b} and the update by \cite{Ghan2011} to \cite{Abdul-Razzak2000}, respectively; the curves correspond to 4th order chaos expansions with coefficients derived using the named method.} 
    \label{fig:single_var_sensitivity}
\end{figure}

<!--
Figure 4: Stratified error plots, looking at how the parameterization  performs for different levels of pollution and updraft speed
-->
\begin{figure}[t]
    \noindent
    \includegraphics[width=\textwidth]{figures/2_SENSITIVITY/updraft_pollution_act_frac_error.pdf}
    \caption{Mean relative error in activated fraction for 4th-order OLS-derived chaos expansion (a), ARG (b), and MBN (c) schemes relative to detailed parcel model. Updraft speeds range from light (10-50 \si{\cm\per\second}), moderate (0.5-2.0 \si{\metre\per\second}), and strong (2.0-10.0 \si{\metre\per\second}); pollution levels range from clean (10-250 \si{\per\cubic\cm}), light (250-1000 \si{\per\cubic\cm}), moderate (1000-2500 \si{\per\cubic\cm}), and heavy (2500-10000 \si{\per\cubic\cm}). Error bars denote 95\% confidence interval on mean relative error from in-bin samples; samples with $\mu < 10\,\si{\nano\metre}$ were omitted from these calculations.}
    \label{fig:stratified_error}
\end{figure}

<!--
![Test. This is the caption! \label{fig:nenes_comp} ](nenes_comp_vs_obs.pdf)

![This is a second test. \label{fig:example_pm}](model_example.png)
-->
