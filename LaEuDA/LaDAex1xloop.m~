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
np = 10; nqq = 10; % np is the total number of particles in the 
beta = 1; % ocean drag coefficients; 1/apha gives the decorrelation time as in a OU-process

domain = [-pi pi -pi pi];
nx = 2 * K_max + 1; ny = nx; ndim = nx^2;

%% get the number density data; 230400=480^2; 129600=360^2
maxo = solveParticleModel(domain, sigma_xy, sigma_v, np, dt, kk, rk, N, u_hat,beta);
toc

%% LaDA data
for nqq = 2:np
    %% LaDA data
    nqq
    x = zeros(nqq,N); y = zeros(nqq,N);
    %np = 25600*16; % number of total particles; number of drifts/sparse  data for DA
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
    % a matrix used in the filtering formulae
    InvBoB = eye(2*nqq)/sigma_xy/sigma_xy;
    n = length(kk(1,:)) + 2*nqq; %
    mu0 = zeros(n, 1); % initial value of posterior mean; u_hat(:,1) is zero
    R0 = zeros(n,n); % initial value of posterior covariance
    u_post_mean = zeros(n,N); % posterior mean
    u_post_mean(:,1) = mu0;
    u_post_cov = zeros(n,N); % posterior covariance
    u_post_cov(:,1) = diag(R0); % only save the diagonal elements
    A1 = zeros(2*nqq,n); A1(:, Dim_U + 1:end) = eye(2*nqq);
    Q = zeros(2*nqq,Dim_U); QQ = zeros(Dim_U, 2*nqq);
    sa0 = [a0; zeros(2*nqq,1)]; % sum a0 part
    b2dotb2 = [Sigma_u*Sigma_u' QQ; Q zeros(2*nqq)];
    for i = 2:N
        % matrix for filtering
        xdiff = [x(:,i)-x(:,i-1); y(:,i)-y(:,i-1)];
        xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
        xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;

        Q(1:nqq,:)      = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:));
        Q(nqq+1:2*nqq,:) = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:));
        Q = beta*Q;
        sa1 = [a1 QQ; Q -beta*eye(2*nqq)];

        % update the posterior mean and posterior covariance
        mu = mu0 + (sa0 + sa1 * mu0) * dt + (R0 * A1') * InvBoB * (xdiff - A1 * mu0 * dt);  % A0 = 0
        R = R0 + (sa1 * R0 + R0* sa1' + b2dotb2 - (R0*A1') * InvBoB * (R0*A1')')*dt;
        u_post_mean(:,i) = mu;
        u_post_cov(:,i) = diag(R);
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


    %% calculate RMSE and PCC in physical domain
    rmsepccPhyDomain
    save(['./err/err' num2str(nqq,'%05.f') '.mat'], "rmse", "rrmse", "pcc")

    save(['./uhat/uhat' num2str(nqq,'%05.f') '.mat'], "u_post_mean")
end


%% The following lines are for plotting the results
figure
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    %indd = mod(24*(i-1)+1,40); % for kmax = 3
    indd = mod(60*(i-1)+1,100); % for kmax = 5
    plot(dt:dt:N*dt, u_hat(indd,1:N), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, u_post_mean(indd,:), 'r', 'linewidth',2)
    title(['(a) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(indd,:))+2*sqrt(real(u_post_cov(indd,:))), real(u_post_mean(indd,end:-1:1))-2*sqrt(real(u_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')


    subplot(4,2,2*i)
    hold on
    indd = Dim_Ug*2 + 15*(i-1)+4;
    plot(dt:dt:N*dt, u_hat(indd,1:N), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, u_post_mean(indd,:), 'r', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(indd,:))+2*sqrt(real(u_post_cov(indd,:))), real(u_post_mean(indd,end:-1:1))-2*sqrt(real(u_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(vx, pny, pnx);
vy = reshape(vy, pny, pnx);
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])


return

%%

figure
for i = 1:4
    subplot(2,2,i)
    if i == 1
        hold on
        plot(dt:dt:N*dt, u_hat(1,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, u_post_mean(1,:), 'r', 'linewidth',2)
        title(['(a) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
        patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(1,:))+2*sqrt(real(u_post_cov(1,:))), real(u_post_mean(1,end:-1:1))-2*sqrt(real(u_post_cov(1,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    elseif i == 2
        hold on
        plot(dt:dt:N*dt, u_hat(Dim_Ug*2+1,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, u_post_mean(Dim_Ug*2+1,:), 'r', 'linewidth',2)
        patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(Dim_Ug*2+1,:))+2*sqrt(real(u_post_cov(Dim_Ug*2+1,:))), real(u_post_mean(Dim_Ug*2+1,end:-1:1))-2*sqrt(real(u_post_cov(Dim_Ug*2+1,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
        title(['(b) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 3
        hold on
        plot(dt:dt:N*dt, u_hat(6,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, u_post_mean(6,:), 'r', 'linewidth',2)
        patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(6,:))+2*sqrt(real(u_post_cov(6,:))), real(u_post_mean(6,end:-1:1))-2*sqrt(real(u_post_cov(6,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
        title(['(c) GB mode ( ', num2str(kk(1,6)),' , ', num2str(kk(2,6)), ' )'],'fontsize',14)
    elseif i == 4
        hold on
        plot(dt:dt:N*dt, u_hat(Dim_Ug*2+6,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, u_post_mean(Dim_Ug*2+6,:), 'r', 'linewidth',2)
        patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(Dim_Ug*2+6,:))+2*sqrt(real(u_post_cov(Dim_Ug*2+6,:))), real(u_post_mean(Dim_Ug*2+6,end:-1:1))-2*sqrt(real(u_post_cov(Dim_Ug*2+6,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
        title(['(d) GB mode ( ', num2str(kk(1,6)),' , ', num2str(kk(2,6)), ' )'],'fontsize',14)
    end
    set(gca,'fontsize',12)
    box on
    xlabel('t')
end

%%
figure
ss = 1;
for i = 2:2:round(T/dt/100)

    u = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([0, 2*pi ])
    ylim([0, 2*pi ])
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    box on
    title(['t = ', num2str(dt*100*(i-1))])
    hold on
    plot(x(:,1+100*(i-1)),y(:,1+100*(i-1)),'ko')
    pause(0.1);
    hold off
    ss = ss + 1;
    set(gca,'fontsize',12)
end

%%
figure
for i = 1:9
    subplot(3,3,i)
    hold on
    u = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([0, 2*pi ])
    ylim([0, 2*pi ])
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    box on
    title(['t = ', num2str(dt*100*(i-1))])
    plot(x(:,1+100*(i-1)),y(:,1+100*(i-1)),'ko','linewidth',2)
    set(gca,'fontsize',12)
end