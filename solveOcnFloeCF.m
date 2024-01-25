% Lagrangian data assimilation with a given number of tracer L with
% observation of velocity
close all; clc; clear all; 
tic
rng(11); % fix the random number seed to reproduce results

% generate ocean current
OU_SWEex3 % kmax = 6; incompressible flow; only GB modes
timeOcn = toc

sigma_xy = 0.001; % noise in the Lagrangian tracer equations
sigma_v = 0.01;
np = 2916;  % np is the total number of particles in the 
beta = 0.1; % ocean drag coefficients; 1/apha gives the decorrelation time as in a OU-process

domain = [-pi pi -pi pi];
nx = 2 * K_max + 1; ny = nx; ndim = nx^2; hx = (domain(2) - domain(1))/nx; hy = hx;

%% get the number density data; 230400=480^2; 129600=360^2
maxo = solveParticleModelCF(domain, sigma_xy, sigma_v, np, dt, kk, rk, N, u_hat,beta);
timeParticle = toc

%save('./data/ws.mat')
save('/g/data/zv32/seaIceFloeData/datakmax4beta01N2196/ws.mat');

