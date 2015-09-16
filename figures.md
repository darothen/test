
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
