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
close all; clc; clear all;
rng(77); % fix the random number seed to reproduce results
tic
beta = 1;
domain = [-pi pi -pi pi];
nx = 9; ny = nx; ndim = nx^2; hx = (domain(2) - domain(1))/nx; hy = (domain(4) - domain(3))/ny;

OU_SWEex3 % incompressible flow; only GB modes
timeEuDAocn = toc
% load('./uhat/ocn.mat', "u_hat","kk","rk");
% 
sigma_xy = 0.001; % noise in the Lagrangian tracer equations
sigma_v = 0.0;

np = 2000;
maxo = solveParticleModelCF(domain, sigma_xy, sigma_v, np, dt, kk, rk, N, u_hat,beta);

%% LaDA data
nqq = 40;
x = zeros(nqq,N); y = zeros(nqq,N); 
aa = randperm(np); aa = aa(1:nqq);
for j=1:N % observational data
    file_name = sprintf('./data/time%05d.mat', j);
    load(file_name)

    x(:,j) = La(aa,1); y(:,j) = La(aa,2); 
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Full filter & Smoother %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Full filter')

setsigv = (1:2)*0.02 - 0.02;

for nsv = 1:length(setsigv)
    
    sigma_v = setsigv(nsv);
    
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
    %b2dotb2 = [zeros(2*nqq) Q; QQ b2b];
    b2dotb2 = [eye(2*nqq)*sigma_v*sigma_v Q; QQ b2b];
    % u_post_cov_all = zeros(n,n,N);
    % u_post_cov_all(:,:,1) = R0;
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
        %u_post_cov_all(:,:,i) = R;
        
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
    
    timeLaDA = toc
    save(['./uhat/nsv' num2str(nsv,'%02.f') 'uhat.mat'],"u_post_mean", "u_post_cov");

    %%
    rmsepccPhyDomain4LaDA

end

%%
errpcc = zeros(2, length(setsigv));
for nsv = 1:length(setsigv)
    
    load(['./err/nsv' num2str(nsv,'%01.f') 'err.mat'], "err")

    errpcc(1, nsv) = mean(0.5*(err(20000:end, 1) + err(20000:end, 2) ) );
    errpcc(2, nsv) = mean(0.5*(err(20000:end, 3) + err(20000:end, 4) ) );
end

%%

figure
subplot(2,3,3)
hold on
yyaxis left
plot(setsigv, errpcc(1,:), '-*', 'linewidth',2)
title('Normalised RMSE','fontsize',24)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on

hold on
yyaxis right
plot(setsigv, errpcc(2,:), '--s', 'linewidth',2)
title('PCC','fontsize',24)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on
xlabel('noise v')




return


%% The following lines are for plotting the results
load(['./uhat/nsv' num2str(6,'%02.f') 'uhat.mat'],"u_post_mean", "u_post_cov");
%figure
subplot(2,3,2)
hold on
%indd = mod(24*(i-1)+1,40); % for kmax = 3
indd = 20;
plot(dt:dt:N*dt, real(u_hat(indd,1:N)), 'b', 'linewidth',2)
plot(dt:dt:N*dt, real(u_post_mean(2*nqq+indd,:)), 'r', 'linewidth',2)
title(['GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',24)
patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq+indd,:))+2*sqrt(real(u_post_cov(2*nqq+indd,:))), real(u_post_mean(2*nqq+indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq+indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
set(gca,'fontsize',24); set(gca,'linewidth',2)

box on
%xlabel('t')


subplot(2,3,5)
hold on
indd = 2;
plot(dt:dt:N*dt, real(u_hat(indd,1:N)), 'b', 'linewidth',2)
plot(dt:dt:N*dt, real(u_post_mean(2*nqq+indd,:)), 'r', 'linewidth',2)
patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq+indd,:))+2*sqrt(real(u_post_cov(2*nqq+indd,:))), real(u_post_mean(2*nqq+indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq+indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
title(['GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',24)
set(gca,'fontsize',24); set(gca,'linewidth',2)
set(gca,'linewidth',2)
box on
xlabel('t')


%%

figure
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    %indd = mod(24*(i-1)+1,40); % for kmax = 3
    indd = mod(60*(i-1)+1,100); % for kmax = 5
    plot(dt:dt:N*dt, real(u_hat(indd,1:N)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean(2*nqq+indd,:)), 'r', 'linewidth',2)
    title(['(a) Gravity mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq+indd,:))+2*sqrt(real(u_post_cov(2*nqq+indd,:))), real(u_post_mean(2*nqq+indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq+indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = Dim_Ug*2 + 15*(i-1)+4;
    plot(dt:dt:N*dt, real(u_hat(indd,1:N)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean(2*nqq+indd,:)), 'r', 'linewidth',2)
    title(['(a) Gravity mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq+indd,:))+2*sqrt(real(u_post_cov(2*nqq+indd,:))), real(u_post_mean(2*nqq+indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq+indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
end

%%
figure
pnx = 32; pny = 32;
[nxx,nyy] = meshgrid(linspace(-pi,pi,pnx), linspace(-pi,pi,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 1000;
subplot(2,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('SWE current')

subplot(2,4,5)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 10000;
subplot(2,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,6)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 30000;
subplot(2,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,7)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])


%--------------------------
ind = 50000;
subplot(2,4,4)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,8)
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LaDA')
xlabel(['t = ', num2str(dt*ind)])

save('./ws/ws.mat')
