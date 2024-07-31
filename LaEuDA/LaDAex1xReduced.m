%********************************************************************
%
%  LEMDA: Lagrangian Eulerian Multiscale Data Assimilation (LaDA Part)
%
%  If you intend to use this package, please acknowledge the authors and the
%  source, particularly by citing our work on LEMDA.
%
%  Permission to use this work is granted, provided that the copies
%  are not made or distributed for commercial purposes and provided that
%  NO WARRANTY, express or implied, of merchantability or fitness for any
%  particular purpose is made with respect to products which have been altered,
%  subjected to abuse or improperly used. Under no circumstances will the
%  author be liable for any damages, lost profits, malpractice claims, 
%  inconvenience, or any other losses caused by the use of this package.
%
%  If anyone needs other permissions that aren't covered by the above,
%  please contact the authors. All rights reserved.
%
%  (c) Quanling Deng, 2022-2024
%  Quanling.Deng@anu.edu.au; qdeng12@gmail.com
%
%  School of Computing
%  Australian National University
%
%********************************************************************

% Lagrangian data assimilation with a given number of tracer L with
% observation of velocity
close all; clc; clear all; 
tic
rng(11); % fix the random number seed to reproduce results

% generate ocean current
OU_SWEex3 % kmax = 6; incompressible flow; only GB modes

sigma_xy = 0.001; % noise in the Lagrangian tracer equations
sigma_v = 0.00;
np = 20; nqq = 10; % np is the total number of particles in the 
beta = 1; % ocean drag coefficients; 1/apha gives the decorrelation time as in a OU-process

domain = [-pi pi -pi pi];
nx = 2 * K_max + 1; ny = nx; ndim = nx^2;

%% get the number density data; 230400=480^2; 129600=360^2
maxo = solveParticleModel(domain, sigma_xy, sigma_v, np, dt, kk, rk, N, u_hat,beta);
toc

%% LaDA data
x = zeros(nqq,N); y = zeros(nqq,N); 
aa = randperm(np); aa = aa(1:nqq);
for j=1:N % observational data
    file_name = sprintf('./data/time%05d.mat', j);
    load(file_name)

    x(:,j) = La(aa,1); y(:,j) = La(aa,2); 
end
toc

l = length(k(1,:)); % number of Fourier wavenumbers
if sigma_g ~=0
    sgm = [sigma_g^2 * ones(1,2*l), sigma_B^2 * ones(1,l-1), sigma_x_b^2, sigma_y_b^2];
    dp = [d_g* ones(1,2*l), d_B* ones(1,l-1), d_b, d_b];
    R_eq  = diag(sgm/2*2./dp);
    mu_eq = zeros(3*l+1,1);
    Dim = 3*l+1;
else
    sgm = [sigma_B^2 * ones(1,l-1), sigma_x_b^2, sigma_y_b^2];
    dp = [d_B* ones(1,l-1), d_b, d_b];
    R_eq  = diag(sgm/2*2./dp);
    mu_eq = zeros(l+1,1);
    Dim = l+1;
end
% quantify the uncertainty reduction using relative entropy
Relative_Entropy_Signal = zeros(1,N);
Relative_Entropy_Dispersion = zeros(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Full filter & Smoother %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Full filter')
% a matrix used in the filtering formulae
InvBoB = eye(2*nqq)/sigma_xy/sigma_xy;
n = length(kk(1,:)) + 2*nqq; %
mu0 = zeros(n, 1); % initial value of posterior mean; u_hat(:,1) is zero
R0 = zeros(n,n); % initial value of posterior covariance
u_post_mean = zeros(n,N); % posterior mean
u_post_mean(:,1) = mu0;
u_post_cov = zeros(n,N); % posterior covariance
u_post_cov(:,1) = diag(R0); % only save the diagonal elements
A1 = [eye(2*nqq), zeros(2*nqq, Dim_U)];
Q = zeros(2*nqq,Dim_U); QQ = zeros(Dim_U, 2*nqq);
sa0 = [zeros(2*nqq,1); a0]; % sum a0 part
b2b = Sigma_u*Sigma_u';
b2dotb2 = [zeros(2*nqq) Q; QQ b2b];

u_post_cov_all = zeros(n,n,N);
u_post_cov_all(:,:,1) = R0;

for i = 2:N
    % matrix for filtering
    xdiff = [x(:,i)-x(:,i-1); y(:,i)-y(:,i-1)]; 
    xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
    xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;
    
    Q(1:nqq,:)      = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:));
    Q(nqq+1:2*nqq,:) = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:));
    Q = beta*Q;
    sa1 = [-beta*eye(2*nqq) Q; QQ a1]; 
    
    % update the posterior mean and posterior covariance
    mu = mu0 + (sa0 + sa1 * mu0) * dt + (R0 * A1') * InvBoB * (xdiff - A1 * mu0 * dt);  % A0 = 0
    R = R0 + (sa1 * R0 + R0* sa1' + b2dotb2 - (R0*A1') * InvBoB * (R0*A1')')*dt;
    u_post_mean(:,i) = mu;
    u_post_cov(:,i) = diag(real(R));
    u_post_cov_all(:,:,i) = R;

    mu0 = mu;
    R0 = R;
    if sigma_g ~=0
        mu_t = mu;
        R_t = R;
    else
        mu_t = mu(2*l+1:end);
        R_t = R(2*l+1:end,2*l+1:end);
    end
    % computing the information gain via relative entropy
%     Relative_Entropy_Signal(i) = real(1/2 * ( (mu_t-mu_eq)' / R_eq * (mu_t-mu_eq) ));
%     Relative_Entropy_Dispersion(i) = real(1/2 * ( trace(R_t/R_eq) - Dim - log(det( R_t/R_eq)) ));
end
% Relative_Entropy_Signal_All = mean(Relative_Entropy_Signal(1000:end));
% Relative_Entropy_Dispersion_All = mean(Relative_Entropy_Dispersion(1000:end));

disp('Full smoother')


R_end = R;
disp('Numerical value of the diagonal entry of R22 at final time')
disp(R(end,end))
R21 = R(2*nqq+1:end, 1:2*nqq); R12 = R21';
temp = diag(diag(R21*R12/sigma_xy^2));

%%
mu_s = zeros(n,N); 
mu_s(:,end) = mu;
R_s = zeros(n,N); % posterior covariance
R_s(:,1) = diag(R); % only save the diagonal elements
R_s_temp0 = R;

for i = N-1:-1:1
    Q(1:nqq,:)     = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:)));
    Q(nqq+1:2*nqq,:) = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:)));
    sa1 = [-beta* eye(2*nqq,2*nqq), Q; zeros(Dim_U,2*nqq), a1];
    
    mu_s(:,i) = mu_s(:,i+1) + (- sa1*mu_s(:,i+1) + b2dotb2 * (squeeze(u_post_cov_all(:,:,i+1)) \ (u_post_mean(:,i+1) - mu_s(:,i+1)))) * dt;
    
    R_s_temp = R_s_temp0 +  ( - ( sa1 + b2dotb2 / (squeeze(u_post_cov_all(:,:,i+1))) ) * R_s_temp0 ...
        - R_s_temp0 * ( sa1' + (squeeze(u_post_cov_all(:,:,i+1)))\(b2dotb2)  ) + b2dotb2) * dt;

    R_s_temp = (R_s_temp + R_s_temp')/2;
    R_s(:,i) = diag(real(R_s_temp));
    R_s_temp0 = R_s_temp;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Reduced filter & smoother %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reduced filter for R')
search_r_equil
r2_rd = real(r2_rd);

mu0 = zeros(n, 1); % initial value of posterior mean
R0 = 0.01*eye(n,n) + 0.001*rand(n,n); % initial value of posterior covariance (this cannot be zero if smoothing is later applied!)
u_post_mean_rd = zeros(n,N); % posterior mean
u_post_mean_rd(:,1) = mu0;
u_post_cov_rd = zeros(n,N); % posterior covariance
u_post_cov_rd(:,1) = diag(R0); % only save the diagonal elements
u_post_cov_all_rd = zeros(n,n,N);
u_post_cov_all_rd(:,:,1) = R0;
R22 = r2_rd*eye(Dim_U); 

for i = 2:N

    xdiff = [x(:,i)-x(:,i-1); y(:,i)-y(:,i-1)]; 
    xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
    xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;

    Q(1:nqq,:)      = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:));
    Q(nqq+1:2*nqq,:) = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:));
    Q = beta*Q;
    sa1 = [-beta*eye(2*nqq) Q; QQ a1]; 
    
    % update the posterior mean and posterior covariance
    R11_0 = R0(1:2*nqq, 1:2*nqq);
    R21_0 = R0(2*nqq+1:end, 1:2*nqq); R12_0 = R21_0';
    R22_0 = R0(2*nqq+1:end, 2*nqq+1:end);
    %R22_temp = (b2dotb2./(-a1 + sqrt(-a1 + L/sigma_xy^2*b2b )));
    %R22 = diag(diag(R22_temp));
    %R22 = r2_rd*eye(Dim_U); %(- temp + b2b)/2 * inv(-a1);

    R11 = R11_0 + (-2*beta * R11_0 + Q * R21_0 + R12_0 * Q' - 1/sigma_xy^2 * R11_0 * R11_0) * dt;
    R12 = R12_0 + (-beta*R12_0 + Q* R22_0- R12_0* (-a1)- 1/sigma_xy^2 * R11_0 * R12_0) * dt;
    R21 = R12';
    R = [R11,R12;R21,R22];
    %R =  R0 + (sa1 * R0 + R0* sa1' + b2dotb2' - (R0*A1') * InvBoB * (R0*A1')')*dt;
    %R(1:2*nqq,1:2*nqq) = R(1:2*nqq,1:2*nqq) .* eye(2*nqq);
    %R(2*nqq+1:end, 2*nqq+1:end) = R(2*nqq+1:end, 2*nqq+1:end).*eye(n);
    mu = mu0 + (sa0 + sa1 * mu0) * dt + (R0 * A1') * InvBoB * (xdiff - A1 * mu0 * dt);  % A0 = 0
 
    u_post_mean_rd(:,i) = mu;
    u_post_cov_rd(:,i) = diag(real(R));
    mu0 = mu;
    R0 = R;
    u_post_cov_all_rd(:,:,i) = R;
end
R_end_rd = R;

mu_s_rd = zeros(n,N); 
mu_s_rd(:,end) = mu;
R_s_rd = zeros(n,N); % posterior covariance
R_s_rd(:,1) = diag(R); % only save the diagonal elements
R_s_temp0 = R;
for i = N-1:-1:1
    
    Q(1:nqq,:)       = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:)));
    Q(nqq+1:2*nqq,:) = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:)));
    sa1 = [-beta* eye(2*nqq,2*nqq), Q; zeros(Dim_U,2*nqq), a1];

    mu_s_rd(:,i) = mu_s_rd(:,i+1) + (- sa1*mu_s_rd(:,i+1) + b2dotb2 * (squeeze(u_post_cov_all_rd(:,:,i+1)) \ (u_post_mean_rd(:,i+1) - mu_s_rd(:,i+1)))) * dt;
    R_s_temp = R_s_temp0 +  ( - ( sa1 + b2dotb2 / (squeeze(u_post_cov_all_rd(:,:,i+1))) ) * R_s_temp0 ...
        - R_s_temp0 * ( sa1' + (squeeze(u_post_cov_all_rd(:,:,i+1)))\(b2dotb2) ) + b2dotb2) * dt;
    R_s_temp = (R_s_temp + R_s_temp')/2;
    R_s_rd(:,i) = diag(real(R_s_temp));
    R_s_temp0 = R_s_temp;
end


%% plots
figure
for i = [1, 3]
    subplot(2,2,i)
    if i == 1
        hold on
        plot(dt:dt:N*dt, real(u_hat(1,:)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean(2*nqq+1,:)), 'r', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean_rd(2*nqq+1,:)), 'g', 'linewidth',2)
        title(['(a) Filtering mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)        
        post_upper = real(u_post_mean(2*nqq+1,:)) + 2 * sqrt(u_post_cov(2*nqq+1,:));
        post_lower = real(u_post_mean(2*nqq+1,:)) - 2 * sqrt(u_post_cov(2*nqq+1,:));
        post_upper_rd = real(u_post_mean_rd(2*nqq+1,:)) + 2 * sqrt(u_post_cov_rd(2*nqq+1,:));
        post_lower_rd = real(u_post_mean_rd(2*nqq+1,:)) - 2 * sqrt(u_post_cov_rd(2*nqq+1,:));
%     elseif i == 2
%         hold on
%         plot(dt:dt:N*dt, real(u_hat(2,:)), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, real(u_post_mean(2*nqq+2,:)), 'r', 'linewidth',2)
%         plot(dt:dt:N*dt, real(u_post_mean_rd(2*nqq+2,:)), 'g', 'linewidth',2)
%         title(['(b) Filtering mode ( ', num2str(kk(1,2)),' , ', num2str(kk(2,2)), ' )'],'fontsize',14)
%         post_upper = real(u_post_mean(2*nqq+2,:)) + 2 * sqrt(u_post_cov(2*nqq+2,:));
%         post_lower = real(u_post_mean(2*nqq+2,:)) - 2 * sqrt(u_post_cov(2*nqq+2,:));
%         post_upper_rd = real(u_post_mean_rd(2*nqq+2,:)) + 2 * sqrt(u_post_cov_rd(2*nqq+2,:));
%         post_lower_rd = real(u_post_mean_rd(2*nqq+2,:)) - 2 * sqrt(u_post_cov_rd(2*nqq+2,:));
    elseif i == 3
        hold on
        plot(dt:dt:N*dt, real(u_hat(1,:)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(mu_s(2*nqq+1,:)), 'r', 'linewidth',2)
        plot(dt:dt:N*dt, real(mu_s_rd(2*nqq+1,:)), 'g', 'linewidth',2)
        title(['(b) Smoothing mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)        
        post_upper = real(mu_s(2*nqq+1,:)) + 2 * sqrt(R_s(2*nqq+1,:));
        post_lower = real(mu_s(2*nqq+1,:)) - 2 * sqrt(R_s(2*nqq+1,:));
        post_upper_rd = real(mu_s_rd(2*nqq+1,:)) + 2 * sqrt(R_s_rd(2*nqq+1,:));
        post_lower_rd = real(mu_s_rd(2*nqq+1,:)) - 2 * sqrt(R_s_rd(2*nqq+1,:));
%     elseif i == 4
%         hold on
%         plot(dt:dt:N*dt, real(u_hat(2,:)), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, real(mu_s(2*nqq+2,:)), 'r', 'linewidth',2)
%         plot(dt:dt:N*dt, real(mu_s_rd(2*nqq+2,:)), 'g', 'linewidth',2)
%         title(['(d) Smoothing mode ( ', num2str(kk(1,2)),' , ', num2str(kk(2,2)), ' )'],'fontsize',14)        
%         post_upper = real(mu_s(2*nqq+2,:)) + 2 * sqrt(R_s(2*nqq+2,:));
%         post_lower = real(mu_s(2*nqq+2,:)) - 2 * sqrt(R_s(2*nqq+2,:));
%         post_upper_rd = real(mu_s_rd(2*nqq+2,:)) + 2 * sqrt(R_s_rd(2*nqq+2,:));
%         post_lower_rd = real(mu_s_rd(2*nqq+2,:)) - 2 * sqrt(R_s_rd(2*nqq+2,:));
    end
    
    tt = dt:dt:N*dt;
    patch([tt,tt(end:-1:1)],[post_lower,post_upper(end:-1:1)],'r','facealpha',0.2,'linestyle','none');
    patch([tt,tt(end:-1:1)],[post_lower_rd,post_upper_rd(end:-1:1)],'g','facealpha',0.2,'linestyle','none');
    set(gca,'fontsize',12)
    if i == 1
        legend('Truth','Full','Reduced','2std Full','2std Reduced')
        yylim = get(gca,'ylim');
    end
    box on
    xlabel('t')
    xlim([1,9])
    ylim([yylim(1),yylim(2)])
end
subplot(2,2,2)
pcolor(real(R_end))
colorbar
title('(c) Full cov at t = 10 (filter)')
box on
set(gca,'fontsize',12)
axis equal
c = get(gca,'clim');
hold on 
plot([2*nqq,2*nqq],[2*nqq,2*nqq+n],'r')
plot([2*nqq,2*nqq+n],[2*nqq,2*nqq],'r')
subplot(2,2,4)
pcolor(real(R_end_rd))
colorbar
title('(d) Reduced-order cov at t = 10 (filter)')
box on
set(gca,'fontsize',12)
axis equal
caxis([c(1),c(2)])
hold on 
plot([2*nqq,2*nqq],[2*nqq,2*nqq+n],'r')
plot([2*nqq,2*nqq+n],[2*nqq,2*nqq],'r')

%%
figure
for i = 1:4
    subplot(2,4,i)    
    hold on
    plot(dt:dt:N*dt, real(x(i,:)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean(i,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean_rd(i,:)), 'g', 'linewidth',2)
    title(['Filtering mode ( ', num2str(kk(1,i)),' , ', num2str(kk(2,i)), ' )'],'fontsize',14)      
    box on
    set(gca,'fontsize',12)
    if i == 1
        legend('Truth','Full','Reduced')
    end
    xlim([1,9])
    subplot(2,4,i+4)    
    hold on
    plot(dt:dt:N*dt, real(x(i,:)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(mu_s(i,:)), 'r--', 'linewidth',2)
    plot(dt:dt:N*dt, real(mu_s_rd(i,:)), 'g--', 'linewidth',2)
    title(['Smoothing mode ( ', num2str(kk(1,i)),' , ', num2str(kk(2,i)), ' )'],'fontsize',14)      
    box on
    set(gca,'fontsize',12)
    xlim([1,9])
end



return
%% calculate RMSE and PCC in physical domain
rmsepccPhyDomain


%% The following lines are for plotting the results
figure
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    %indd = mod(24*(i-1)+1,40); % for kmax = 3
    indd = mod(60*(i-1)+1,100); % for kmax = 5
    plot(dt:dt:N*dt, u_hat(indd,1:N), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, u_post_mean(2*nqq + indd,:), 'r', 'linewidth',2)
    title(['(a) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq + indd,:))+2*sqrt(real(u_post_cov(2*nqq + indd,:))), real(u_post_mean(2*nqq + indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq + indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = Dim_Ug*2 + 15*(i-1)+4;
    plot(dt:dt:N*dt, u_hat(indd,1:N), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, u_post_mean(2*nqq + indd,:), 'r', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq + indd,:))+2*sqrt(real(u_post_cov(2*nqq + indd,:))), real(u_post_mean(2*nqq + indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq + indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end

figure
pnx = 40; pny = 40;
[nxx,nyy] = meshgrid(linspace(-pi,pi,pnx), linspace(-pi,pi,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 100;
subplot(2,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('SWE current')

subplot(2,4,5)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 500;
subplot(2,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,6)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 1000;
subplot(2,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,7)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])


%--------------------------
ind = 2000;
subplot(2,4,4)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,8)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq + 1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])


