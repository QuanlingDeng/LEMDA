# LEMDA
LEMDA: Lagrangian-Eulerian Multiscale Data Assimilation

This MATLAB package documents the algorithms for the LEMDA framework. In particular, it generates the Figures in the submitted LDEMA paper. Several comments are

- For Figure 2, run LaDAex1xloop.m in folder LaEuDA. There are several loops for different beta values (drag coefficients). We modify the beta values in line 39 for different values used for Panels (a-c) in Figure 1.  Other panels in Figure 1 are special cases of these loops. 

- For Figure 3, run LaDAex1xReduced.m in folder LaEuDA.

- For Figure 5 and panels (a)-(c) of Figure 6, run EuDAloop2.m in folder LaEuDA to generate all the synthetic data and EuDA posterior means and covariances for L=j*500, j=1,2,...16, number of floes. This will produce panel (c) of Figure 6. Then uncomment lines 202 to 373 in EuDAloop2.m to generate the panels (a)-(c) of Figure 5. For panel (d) of Figure 5, run rmsepccPhyDomainEuDA.m and then ``Run Section" (Matlab functionality) of the section of lines 43 to 67. Similar to panel (d) of Figure 5, run the saved data to generate panels (a)-(b) of Figure 6.

- For panel (d) of Figure 6, it is similar to the loop on the number of particles. Herein, the loop is on the grid size. Run EuDAloop3nx.m in folder LaEuDA.

- For Figures 8 and 9, the simulation settings are very similar to those of Figures 5 and 6, except that herein particle collisions are included, leading to a much higher demand on computational hours. To generate these Figures, run EuDAloop2CF.m and EuDAloop3nxCF.m in folder LaEuDA. The simulations for Figures 8 and 9 in the paper were done in the Australian National Computing Infrastructure (NCI) national facility under project zv32 in a parallel environment. There were commented lines in EuDAloop2CF.m and EuDAloop3nxCF.m for the parallel computation setting which can be easily adjusted in a different parallel environment. As these loop files require substantial computational cost, we suggest adjusting these codes accordingly using parallel computing.

- For Figure 10, run LEMDA.m in the LEMDA folder to generate a set of Figures then use them as subfigures for Figure 9.

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
