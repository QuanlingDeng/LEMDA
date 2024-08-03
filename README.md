# LEMDA
LEMDA: Lagrangian-Eulerian Multiscale Data Assimilation

This MATLAB package documents the algorithms for the LEMDA framework. In particular, it generates the Figures in the submitted LDEMA paper. Several comments are

- For Figure 9, run LEMDA.m in the LEMDA folder to generate a set of Figures then use them as subfigures for Figure 9.

- For Figure 8, run EuDAloop2CF.m and EuDAloop3nxCF.m in folder LaEuDA. These loop files require substantial computational cost, there were run in clusters. We suggest using parallel computing here as well.
  
- For Figure 7, run EuDACF.m in folder LaEuDA; 

- For Figure 6, run LaDAloopSigv.m in folder LaEuDA

- For Figure 5, run EuDAloop2CF.m and EuDAloop3nxCF.m in folder LaEuDA.

- For Figure 4, run EuDAloop2.m and EuDAloop3nx.m in folder LaEuDA.

- For Figure 3, run LaDAex1xReduced.m in folder LaEuDA.

- For Figure 2, run LaDAex1xloop.m in folder LaEuDA. There are several loops for different beta values (drag coefficients). We modify the beta values in line 39 for different values used for Panels (a-c) in Figure 1.  Other panels in Figure 1 are special cases of these loops. 

********************************************************************

 LEMDA: Lagrangian Eulerian Multiscale Data Assimilation 

 If you intend to use this package, please acknowledge the authors and the
 source, particularly by citing our work on LEMDA.

 Permission to use this work is granted, provided that the copies
 are not made or distributed for commercial purposes and provided that
 NO WARRANTY, express or implied, of merchantability or fitness for any
 particular purpose is made with respect to products which have been altered,
 subjected to abuse or improperly used. Under no circumstances will the
 author be liable for any damages, lost profits, malpractice claims, 
 inconvenience, or any other losses caused by the use of this package.

 If anyone needs other permissions that aren't covered by the above,
 please contact the authors. All rights reserved.

 (c) Quanling Deng, 2022-2024
 Quanling.Deng@anu.edu.au; qdeng12@gmail.com

 School of Computing
 Australian National University

********************************************************************
