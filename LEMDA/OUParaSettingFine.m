% Using linear stochastic models to simulate 2D rotating shallow water 
% equation with random amplitudes
% setting for coarse scale modes
fK_max = 1; % the range of Fourier modes is [-K_max, K_max]^2
fk = zeros(2, (2 * fK_max + 1) * (2 * fK_max + 1)); % Total number of Fourier wavenumbers

% arranging Fourier wavenumbers
% arranging in such a way that the complex conjugates modes are next to
% each other, namely (-k1,-k2) will be next to (k1,k2). This will
% facilitate data assimilation and allows it to be a nearly block diagonal
% matrix where each block is 2 by 2.
fm = 1;
for i = - fK_max : fK_max
    if i < 0
        for j = - fK_max: i
            fk(1, fm) = i;
            fk(2, fm) = j;
            fm = fm + 2;
        end
    else
        for j = - fK_max : i - 1
            fk(1, fm) = i;
            fk(2, fm) = j;
            fm = fm + 2;
        end
    end
end
fk(:, 2: 2: end - 1) = - fk(:, 1 : 2 : end - 2);


%%
% Total number of Fourier modes is 3 times as many as that of the Fourier
% wavenumbers because each Fourier wavenumber has three modes: one GB mode
% that is incompressible and two gravity modes that are compressible
fkk_end = [0;0];
fkk = [fk,fk,fkk_end,fk,fkk_end]; % total number of the waves, 2k gravity, k-1 GB and two background modes
fkk(:,length(fk(1,:))) = []; % putting the two +/-(0,0) modes of the gravity waves together
epsilon = 1; 1;0.2; % Rossby number; 0.2 for fast changing flows (need small dt for resolution) while 1 for slowly changing flows
% dispersion relationship for the two gravity modes, p and m stand for plus and minus
% the dispersion relationship for the GB mode is omega = 0
fomegak_p = 1/epsilon * sqrt(fk(1,:).^2 + fk(2,:).^2 + 1); 
fomegak_m = - 1/epsilon * sqrt(fk(1,:).^2 + fk(2,:).^2 + 1);
fomegak = reshape([fomegak_p; fomegak_m],1,[]);
% % omegak = [omegak_p(1:24), omegak_m(1:24), omegak_p(end), omegak_m(end)];
% eigenvectors of GB and gravity modes
% the last column of rk2 and rk3 are the (0,0) gravity modes, which need to
% deal with in a different way
% rk1: GB; rk2: gravity +; rk3: gravity -.
frk1 = [1./sqrt(fk(1,:).^2 + fk(2,:).^2+1) .* (-1i * fk(2,:));
    1./sqrt(fk(1,:).^2 + fk(2,:).^2+1) .* (1i * fk(1,:))];
frk1(:,end) = [1;0]; rk1_end = [0;1]; frk1 = [frk1,rk1_end];
frk2 = [1./sqrt(fk(1,:).^2 + fk(2,:).^2)/sqrt(2)./sqrt(fk(1,:).^2 + fk(2,:).^2 + 1) .* (1i * fk(2,:) + fk(1,:) .* sqrt(fk(1,:).^2 + fk(2,:).^2 + 1));
    1./sqrt(fk(1,:).^2 + fk(2,:).^2)/sqrt(2)./sqrt(fk(1,:).^2 + fk(2,:).^2 + 1) .* (-1i * fk(1,:) + fk(2,:) .* sqrt(fk(1,:).^2 + fk(2,:).^2 + 1))];
frk2(:,end) = [1i;1]/sqrt(2);
frk3 = -[1./sqrt(fk(1,:).^2 + fk(2,:).^2)/sqrt(2)./sqrt(fk(1,:).^2 + fk(2,:).^2 + 1) .* (1i * fk(2,:) - fk(1,:) .* sqrt(fk(1,:).^2 + fk(2,:).^2 + 1));
    1./sqrt(fk(1,:).^2 + fk(2,:).^2)/sqrt(2)./sqrt(fk(1,:).^2 + fk(2,:).^2 + 1) .* (-1i * fk(1,:) - fk(2,:) .* sqrt(fk(1,:).^2 + fk(2,:).^2 + 1))];
frk3(:,end) = [-1i;1]/sqrt(2);

frk = zeros(size(fkk));
frk(:, 1:2:length(fk(1,:))*2-1) = frk2;
frk(:, 2:2:length(fk(1,:))*2) = frk3;
frk(:, length(fk(1,:))*2+1:end) = frk1;


% stochastic systems for each Fourier mode
% written in the vector form
% dU = (a1 U + a0 ) dt + b1 dW

rng(5000) % fix the random number seed
fDim_U = length(fkk(1,:)); % dimension of the system
fDim_Ug = length(fk(1,:)); fDim_UB = fDim_Ug - 1;
fd_B = 1; %0.5; % damping of the GB modes
fd_g = 0.5; % damping of the gravity modes
fd_b = 0.5; % dampling of the two background modes
fsigma_B = 0.01; % noise of the GB modes
fsigma_g = 0.01; % noise of the gravity modes
f_x_b = 0; % large-scale forcing background in x direction
f_y_b = 0; % large-scale forcing background in y direction
fsigma_x_b = 0.1; % noise of large-scale background flow in x direction
fsigma_y_b = 0.1; % noise of large-scale background flow in y direction

% b1: noise coefficient; a1: damping and phase 
fb1 = zeros(fDim_U, fDim_U);
fa1 = - diag([fd_g * ones(1,fDim_Ug * 2), fd_B * ones(1,fDim_UB), fd_b, fd_b]) + 1i * diag([fomegak, zeros(1,fDim_UB+2)]);
fL_u = fa1;
for i = 1: 2: 2 * fDim_Ug + fDim_UB
    if i <= 2 * fDim_Ug
        fb1(i,i) = 1 / sqrt(2) * fsigma_g;
        fb1(i+1, i+1) = -1i / sqrt(2) * fsigma_g;
        fb1(i, i+1) = 1i / sqrt(2) * fsigma_g;
        fb1(i+1, i) = 1 / sqrt(2) * fsigma_g;
    else
        fb1(i,i) = 1 / sqrt(2) * fsigma_B;
        fb1(i+1, i+1) = -1i / sqrt(2) * fsigma_B;
        fb1(i, i+1) = 1i / sqrt(2) * fsigma_B;
        fb1(i+1, i) = 1 / sqrt(2) * fsigma_B;
    end
end
fSigma_u = fb1;
fSigma_u(end-1,end-1) = fsigma_x_b;
fSigma_u(end,end) = fsigma_y_b;

% numerical integration
frd = zeros(fDim_U,N);
% the following two lines can be commented out if the system has no gravity
% modes
frd(1:2:end-3,:) = randn(fDim_Ug + fDim_UB/2, N) + 1i * randn(fDim_Ug + fDim_UB/2, N);
frd(2:2:end-2,:) = conj(frd(1:2:end-3,:)); % noise needs to be complex conjudate within each 2x2 block
frd(1:end-2, :) = randn(fDim_Ug * 2 + fDim_UB, N); 
fa0 = zeros(fDim_U, 1);


%% time series post-processing to get only GB modes = incompressible flows
fDim_Ug = 0; fDim_U = fDim_UB + 2;
fa0 = fa0(length(fk(1,:))*2+1:end);
fa1 = fa1(length(fk(1,:))*2+1:end, length(fk(1,:))*2+1:end); % a1 is actually a symmetric matrix
fkk = fkk(:, length(fk(1,:))*2+1:end);
frk = frk(:, length(fk(1,:))*2+1:end);
fL_u = fL_u(length(fk(1,:))*2+1:end, length(fk(1,:))*2+1:end); 
fSigma_u = fSigma_u(length(fk(1,:))*2+1:end, length(fk(1,:))*2+1:end); 

