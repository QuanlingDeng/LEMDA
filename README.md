# LEMDA
LEMDA: Lagrangian-Eulerian Multiscale Data Assimilation

This MATLAB package documents the algorithms for the LEMDA framework. In particular, it generates the Figures in the LDEMA paper. Several comments are

- Run LEMDA.m in LEMDA folder to generate a set of Figures which are used as subfigures for Figure 9.

- For Figure 7, run EuDACF.m in folder LaEuDA; for Figure 8, run EuDAloop2CF.m and EuDAloop3nxCF.m in folder LaEuDA. These loop files requires substantial computational cost, there were ran in clusters. We suggest using parallel computing here as well.

- For Figure 6, run LaDAloopSigv.m in folder LaEuDA

- For Figure 5, run EuDAloop2.m and EuDAloop3nx.m in folder LaEuDA.

- For Figure 4, run EuDAloop2.m and EuDAloop3nx.m in folder LaEuDA.

- For Figure 3, run LaDAex1xReduced.m in folder LaEuDA.

- For Figure 2, run LaDAex1xloop.m in folder LaEuDA.
