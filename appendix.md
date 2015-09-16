\appendix

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
        