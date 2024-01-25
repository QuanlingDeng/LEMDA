% Eulerian data assimilation with observational data on number density
%close all; clc; clear all;
tic; rng(101);
poolobj = gcp('nocreate'); % If no pool, do not create new one.
c = parcluster('local');
c.NumWorkers = 16;
%if ~isempty(poolobj); delete(poolobj); end
if isempty(poolobj); parpool(c,c.NumWorkers); parpool.IdleTimeout = 1200000000; end


beta = 1;
domain = [-pi pi -pi pi];

% generate ocean current
OU_SWEex2 % incompressible flow; only GB modes
timeEuDAocn = toc
save('/g/data/zv32/seaIceFloeData/cfdata/ocn.mat', "u_hat","kk","rk");
% file_name = sprintf('./data90kkmax4dt4/uhat/ocn.mat');
% load(file_name)
nx = 9; 1*(2 * K_max + 1); ny = nx; ndim = nx^2;

%% get the number density data; 230400=480^2; 129600=360^2
sigma_xy = 0.001; % noise in the Lagrangian tracer equations
sigma_v = 0.00;

npset=(1:16)*500;
parfor npi=1:length(npset)
    
    np = npset(npi); % nqq = 80; % np is the total number of particles in the; nqq observed
    [~] = solveParticleModelCF(domain, sigma_xy, sigma_v, np, npi, dt, kk, rk, N, u_hat,beta);
end

