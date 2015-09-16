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
