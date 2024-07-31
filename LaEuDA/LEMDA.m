% notes: case 1, ro=1, demonstrate euda better than lada; advantage of EuDA
% (using 10x10) over LaDA; then advantages when using LEMDA
% demonstrate advantages of lada and euda; (lada advantage is clear; computing costs)
% how to choose noise strength; euda advantage 
%close all; clc; clear all; 
% Eulerian data assimilation with observational data on number density
% rng(77); % fix the random number seed to reproduce results
% tic 
% beta = 0.1; 
% domain = [-pi pi -pi pi];
% % 
% OU_SWEex2 % kmax = 6; incompressible flow; only GB modes
% save('./uhat/ocn.mat', "u_hat","kk","rk");
% % timeLaDAocn = toc;
% nx = 2 * K_max + 1; ny = nx; ndim = nx^2;
% 
% % get the number density data; 230400=480^2; 129600=360^2
% sigma_xy = 0.001; % noise in the Lagrangian tracer equations
% sigma_v = 0.01;
% np = 49; % np is the total number of particles in the; nqq observed
% 
% maxo = solveParticleModelCF(domain, sigma_xy, sigma_v, np, dt, kk, rk, N, u_hat,beta);
% timeLEMDAparticle = toc;

%
% nqq = 20;
% x = zeros(nqq,N); y = zeros(nqq,N);
% hx = (domain(2) - domain(1))/nx;
% hy = (domain(4) - domain(3))/ny;
% npsol = zeros(ny, nx, N);
% npusol = zeros(ny, nx, N);
% npvsol = zeros(ny, nx, N);
% randpm = randperm(np);
% aa = randpm(1:nqq); 
% for j=1:N
%     file_name = sprintf('./data/time%05d.mat', j);
%     load(file_name)
%     
%     indx = ceil( ( La(:,1) - domain(1) )/hx );
%     indy = ceil( ( La(:,2)  - domain(3) )/hy );
% 
%     for k=1:np
%         npsol(indy(k), indx(k), j) = npsol(indy(k), indx(k), j) + 1;
%         npusol(indy(k), indx(k), j) = npusol(indy(k), indx(k), j) + La(k,3);
%         npvsol(indy(k), indx(k), j) = npvsol(indy(k), indx(k), j) + La(k,4);
%     end
% 
%     x(:,j) = La(aa,1); y(:,j) = La(aa,2);
% end
% for j=1:N
%     [npsol(:,:,j)] = ProcEuNumData(ny, nx, npsol(:,:,j) );
% end
% npusol = npusol*(2*pi)^2/(np*hx*hy);
% npvsol = npvsol*(2*pi)^2/(np*hx*hy);
% npsol = npsol*(2*pi)^2/(np*hx*hy); % number density scaling to around 1.
% npusol0 = npusol; npvsol0 = npvsol; npsol0 = npsol; 
% % ws = 1000;
% % npusol = smoothdata(npusol0,3,"movmean",ws);
% % npvsol = smoothdata(npvsol0,3,"movmean",ws);
% % npsol = smoothdata(npsol0,3,"movmean",ws);
% toc

% mesh
[nxx,nyy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx)); 
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix


%% data assimilation
% sig = getNoise4EuDA(domain, np, 2*nx, nx, dt, 1)
% sigma_n = sqrt(0.5*sig(1)^2 + 0.5*sig(2)^2) % noise in the number density discretized ODEs
toc

% a matrix used in the filtering formulae
sigma_n = 0.02;
InvBoB = diag( [ones(2*nqq,1)/sigma_xy/sigma_xy;  ones(2*ndim,1)/sigma_n/sigma_n] );
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
sigma_v = 2*sigma_v;
b2dotb2 = [eye(2*nqq)*sigma_v*sigma_v Q; QQ b2b];
for i = 2:N
    
%     sigma_n = 0.002 + 0.02*i*dt;
%     InvBoB = diag( [ones(2*nqq,1)/sigma_xy/sigma_xy;  ones(2*ndim,1)/sigma_n/sigma_n] );

    xdiff = [x(:,i)-x(:,i-1); y(:,i)-y(:,i-1)]; 
    xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
    xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;
    
    Q(1:nqq,:)      = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:));
    Q(nqq+1:2*nqq,:) = exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:));
    Q = beta*Q;
    sa1 = [-beta*eye(2*nqq) Q; QQ a1]; 

    nnsol = npsol(:,:,i-1);
    npsols = reshape(npsol(:,:,i-1),[],1);
    npusols = reshape(npusol(:,:,i-1),[],1);
    npvsols = reshape(npvsol(:,:,i-1),[],1);
    dnn = zeros(2*ndim,1);
    dnn(1:ndim,1) = reshape(npusol(:,:,i) - npusol(:,:,i-1),[],1);   
    dnn(1+ndim:2*ndim,1) = reshape(npvsol(:,:,i) - npvsol(:,:,i-1),[],1);   

    A00 = zeros(2*ndim,1);
    dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1), npvsol(:,:,i-1), npusol(:,:,i-1)./nnsol);
    A00(1:ndim,1) = reshape(dvv,[],1);
    dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1), npvsol(:,:,i-1), npvsol(:,:,i-1)./nnsol);
    A00(1+ndim:2*ndim,1) = reshape(dvv,[],1);
    A00 = A00 + beta*[npusols; npvsols];
    A00 = -A00;

    A2 = zeros(2*ndim,Dim_U);
    A2(1:ndim,:) = beta*npsols.*exp(1i * nx_vec * kk) .* (ones(ndim,1) * rk(1,:));
    A2(1+ndim:2*ndim,:) = beta*npsols.*exp(1i * nx_vec * kk) .* (ones(ndim,1) * rk(2,:));
    
    A3 = [zeros(2*nqq,1); A00]; 
    A = [A1; zeros(2*ndim,2*nqq) A2];
    dn = [xdiff; dnn];

    % update the posterior mean and posterior covariance
    mu = mu0 + (sa0 + sa1 * mu0) * dt + (R0 * A') * InvBoB * (dn - (A3 + A * mu0) * dt);  
    R = R0 + (sa1 * R0 + R0* sa1' + b2dotb2 - (R0*A') * InvBoB * (R0*A')')*dt;
    u_post_mean(:,i) = mu;
    u_post_cov(:,i) = diag(R);
    mu0 = mu;
    R0 = R;  
end
toc

%% The following lines are for plotting the results
%u_post_mean = abs(u_post_mean);
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
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(2*nqq+indd,:))+2*sqrt(real(u_post_cov(2*nqq+indd,:))), real(u_post_mean(2*nqq+indd,end:-1:1))-2*sqrt(real(u_post_cov(2*nqq+indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end

%%
rmsepccPhyDomain2

%%
figure
pnx = 32; pny = 32;
[nxx,nyy] = meshgrid(linspace(-pi,pi,pnx), linspace(-pi,pi,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 1000;
subplot(2,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 2000;
subplot(2,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 5000;
subplot(2,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])


%--------------------------
ind = 9900;
subplot(2,4,4)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
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
vx = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(nxx,nyy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

save('./ws/LEMDA.mat')

