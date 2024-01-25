% Using linear stochastic models to simulate 2D rotating shallow water 
% equation with random amplitudes
% setting for coarse scale modes
cK_max = 4; % the range of Fourier modes is [-K_max, K_max]^2
ck = zeros(2, (2 * cK_max + 1) * (2 * cK_max + 1)); % Total number of Fourier wavenumbers

% arranging Fourier wavenumbers
% arranging in such a way that the complex conjugates modes are next to
% each other, namely (-k1,-k2) will be next to (k1,k2). This will
% facilitate data assimilation and allows it to be a nearly block diagonal
% matrix where each block is 2 by 2.
cm = 1;
for i = - cK_max : cK_max
    if i < 0
        for j = - cK_max: i
            ck(1, cm) = i;
            ck(2, cm) = j;
            cm = cm + 2;
        end
    else
        for j = - cK_max : i - 1
            ck(1, cm) = i;
            ck(2, cm) = j;
            cm = cm + 2;
        end
    end
end
ck(:, 2: 2: end - 1) = - ck(:, 1 : 2 : end - 2);


%%
% Total number of Fourier modes is 3 times as many as that of the Fourier
% wavenumbers because each Fourier wavenumber has three modes: one GB mode
% that is incompressible and two gravity modes that are compressible
ckk_end = [0;0];
ckk = [ck,ck,ckk_end,ck,ckk_end]; % total number of the waves, 2k gravity, k-1 GB and two background modes
ckk(:,length(ck(1,:))) = []; % putting the two +/-(0,0) modes of the gravity waves together
epsilon = 1; 1;0.2; % Rossby number; 0.2 for fast changing flows (need small dt for resolution) while 1 for slowly changing flows
% dispersion relationship for the two gravity modes, p and m stand for plus and minus
% the dispersion relationship for the GB mode is omega = 0
comegak_p = 1/epsilon * sqrt(ck(1,:).^2 + ck(2,:).^2 + 1); 
comegak_m = - 1/epsilon * sqrt(ck(1,:).^2 + ck(2,:).^2 + 1);
comegak = reshape([comegak_p; comegak_m],1,[]);
% % omegak = [omegak_p(1:24), omegak_m(1:24), omegak_p(end), omegak_m(end)];
% eigenvectors of GB and gravity modes
% the last column of rk2 and rk3 are the (0,0) gravity modes, which need to
% deal with in a different way
% rk1: GB; rk2: gravity +; rk3: gravity -.
crk1 = [1./sqrt(ck(1,:).^2 + ck(2,:).^2+1) .* (-1i * ck(2,:));
    1./sqrt(ck(1,:).^2 + ck(2,:).^2+1) .* (1i * ck(1,:))];
crk1(:,end) = [1;0]; rk1_end = [0;1]; crk1 = [crk1,rk1_end];
crk2 = [1./sqrt(ck(1,:).^2 + ck(2,:).^2)/sqrt(2)./sqrt(ck(1,:).^2 + ck(2,:).^2 + 1) .* (1i * ck(2,:) + ck(1,:) .* sqrt(ck(1,:).^2 + ck(2,:).^2 + 1));
    1./sqrt(ck(1,:).^2 + ck(2,:).^2)/sqrt(2)./sqrt(ck(1,:).^2 + ck(2,:).^2 + 1) .* (-1i * ck(1,:) + ck(2,:) .* sqrt(ck(1,:).^2 + ck(2,:).^2 + 1))];
crk2(:,end) = [1i;1]/sqrt(2);
crk3 = -[1./sqrt(ck(1,:).^2 + ck(2,:).^2)/sqrt(2)./sqrt(ck(1,:).^2 + ck(2,:).^2 + 1) .* (1i * ck(2,:) - ck(1,:) .* sqrt(ck(1,:).^2 + ck(2,:).^2 + 1));
    1./sqrt(ck(1,:).^2 + ck(2,:).^2)/sqrt(2)./sqrt(ck(1,:).^2 + ck(2,:).^2 + 1) .* (-1i * ck(1,:) - ck(2,:) .* sqrt(ck(1,:).^2 + ck(2,:).^2 + 1))];
crk3(:,end) = [-1i;1]/sqrt(2);

crk = zeros(size(ckk));
crk(:, 1:2:length(ck(1,:))*2-1) = crk2;
crk(:, 2:2:length(ck(1,:))*2) = crk3;
crk(:, length(ck(1,:))*2+1:end) = crk1;


% stochastic systems for each Fourier mode
% written in the vector form
% dU = (a1 U + a0 ) dt + b1 dW

rng(5000) % fix the random number seed
cDim_U = length(ckk(1,:)); % dimension of the system
cDim_Ug = length(ck(1,:)); cDim_UB = cDim_Ug - 1;
cd_B = 1; %0.5; % damping of the GB modes
cd_g = 0.5; % damping of the gravity modes
cd_b = 0.5; % dampling of the two background modes
csigma_B = 0.02; % noise of the GB modes
csigma_g = 0.02; % noise of the gravity modes
f_x_b = 0; % large-scale forcing background in x direction
f_y_b = 0; % large-scale forcing background in y direction
csigma_x_b = 0.1; % noise of large-scale background flow in x direction
csigma_y_b = 0.1; % noise of large-scale background flow in y direction

% b1: noise coefficient; a1: damping and phase 
cb1 = zeros(cDim_U, cDim_U);
ca1 = - diag([cd_g * ones(1,cDim_Ug * 2), cd_B * ones(1,cDim_UB), cd_b, cd_b]) + 1i * diag([comegak, zeros(1,cDim_UB+2)]);
cL_u = ca1;
for i = 1: 2: 2 * cDim_Ug + cDim_UB
    if i <= 2 * cDim_Ug
        cb1(i,i) = 1 / sqrt(2) * csigma_g;
        cb1(i+1, i+1) = -1i / sqrt(2) * csigma_g;
        cb1(i, i+1) = 1i / sqrt(2) * csigma_g;
        cb1(i+1, i) = 1 / sqrt(2) * csigma_g;
    else
        cb1(i,i) = 1 / sqrt(2) * csigma_B;
        cb1(i+1, i+1) = -1i / sqrt(2) * csigma_B;
        cb1(i, i+1) = 1i / sqrt(2) * csigma_B;
        cb1(i+1, i) = 1 / sqrt(2) * csigma_B;
    end
end
cSigma_u = cb1;
cSigma_u(end-1,end-1) = csigma_x_b;
cSigma_u(end,end) = csigma_y_b;

% numerical integration
crd = zeros(cDim_U,N);
% the following two lines can be commented out if the system has no gravity
% modes
crd(1:2:end-3,:) = randn(cDim_Ug + cDim_UB/2, N) + 1i * randn(cDim_Ug + cDim_UB/2, N);
crd(2:2:end-2,:) = conj(crd(1:2:end-3,:)); % noise needs to be complex conjudate within each 2x2 block
crd(1:end-2, :) = randn(cDim_Ug * 2 + cDim_UB, N); 
ca0 = zeros(cDim_U, 1);


%% time series post-processing to get only GB modes = incompressible flows
cDim_Ug = 0; cDim_U = cDim_UB + 2;
ca0 = ca0(length(ck(1,:))*2+1:end);
ca1 = ca1(length(ck(1,:))*2+1:end, length(ck(1,:))*2+1:end); % a1 is actually a symmetric matrix
ckk = ckk(:, length(ck(1,:))*2+1:end);
crk = crk(:, length(ck(1,:))*2+1:end);
cL_u = cL_u(length(ck(1,:))*2+1:end, length(ck(1,:))*2+1:end); 
cSigma_u = cSigma_u(length(ck(1,:))*2+1:end, length(ck(1,:))*2+1:end); 

