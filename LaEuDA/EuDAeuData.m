% Eulerian data assimilation with observational data on number density
rng(77); % fix the random number seed to reproduce results
tic 
beta = 1; 
domain = [-pi pi -pi pi];

% % generate ocean current
% OU_SWEex2 % incompressible flow; only GB modes
% save('./uhat/ocn.mat', "u_hat","kk","rk");

nx = 2 * K_max + 1; ny = nx; ndim = nx^2;

hx = (domain(2) - domain(1))/nx;
hy = (domain(4) - domain(3))/ny;
nsol = zeros(ny, nx, N);
npsol = zeros(ny, nx, N);
nusol = zeros(ny, nx, N);
nvsol = zeros(ny, nx, N);
npusol = zeros(ny, nx, N);
npvsol = zeros(ny, nx, N);

[xx,yy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx));
sol = 5 + 0.3*cos(0.5*xx).*cos(0.5*yy);
usol = -0.15*cos(0.5*xx).*sin(0.5*yy);
vsol = 0.15*sin(0.5*xx).*cos(0.5*yy);
usol = usol.*sol; vsol = vsol.*sol; 
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix
nsol(:,:,1) = sol;
nusol(:,:,1) = real(usol);
nvsol(:,:,1) = real(vsol);

for i = 2:N % generating the number density data
    uox = exp(1i * nx_vec * kk) * (u_hat(:,i-1) .* transpose(rk(1,:)));
    uoy = exp(1i * nx_vec * kk) * (u_hat(:,i-1) .* transpose(rk(2,:)));
    uox = reshape(real(uox), ny, nx);
    uoy = reshape(real(uoy), ny, nx);
    vx = usol./sol; vy = vsol./sol;
    
    forcex = beta*(uox.*nsol(:,:,i-1) - nusol(:,:,i-1));
    forcey = beta*(uoy.*nsol(:,:,i-1) - nvsol(:,:,i-1));

    ssol = solveNumDensity(domain, nx, ny, dt, vx, vy, sol*0.0, sol, 0.001);
    nsol(:,:,i) = ssol; sol = ssol;
    
    ussol = solveNumDensity(domain, nx, ny, dt, vx, vy, forcex, usol, 0.001);
    nusol(:,:,i) = ussol; usol = ussol;
    
    vssol = solveNumDensity(domain, nx, ny, dt, vx, vy, forcey, vsol, 0.001);
    nvsol(:,:,i) = vssol; vsol = vssol;
end

%% data assimilation
npsol = nsol; npusol = nusol; npvsol = nvsol;
% ws = 3;
% npusol = smoothdata(nsol,3,"movmean",ws);
% npvsol = smoothdata(nusol,3,"movmean",ws);
% npsol = smoothdata(nvsol,3,"movmean",ws);

sigma_n =  0.02;

l = length(kk(1,:)); % number of Fourier wavenumbers
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
InvBoB = eye(2*ndim)/sigma_n/sigma_n;
mu0 = u_hat(:,1); % initial value of posterior mean
n = length(kk(1,:));
R0 = zeros(n,n); % initial value of posterior covariance
u_post_mean = zeros(n,N); % posterior mean
u_post_mean(:,1) = mu0;
u_post_cov = zeros(n,N); % posterior covariance
u_post_cov(:,1) = diag(R0); % only save the diagonal elements
for i = 2:N
%     sigma_n = 0.002 + 0.01*i*dt;
%     InvBoB = eye(2*ndim)/sigma_n/sigma_n;

    % matrix for filtering; A2 = zeros(ndim, n)
    nnsol = npsol(:,:,i-1);
    npsols = reshape(npsol(:,:,i-1),[],1);
    npusols = reshape(npusol(:,:,i-1),[],1);
    npvsols = reshape(npvsol(:,:,i-1),[],1);
    dnn = zeros(2*ndim,1);
    dnn(1:ndim,1) = reshape(npusol(:,:,i) - npusol(:,:,i-1),[],1);   
    dnn(1+ndim:2*ndim,1) = reshape(npvsol(:,:,i) - npvsol(:,:,i-1),[],1); 

    A00 = zeros(2*ndim,1);
    %dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1), npvsol(:,:,i-1), npusol(:,:,i-1)./nnsol);
    dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1)./nnsol, npvsol(:,:,i-1)./nnsol, npusol(:,:,i-1));
    dvu = dvv;
    A00(1:ndim,1) = reshape(dvv,[],1);
    %dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1), npvsol(:,:,i-1), npvsol(:,:,i-1)./nnsol);
    dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1)./nnsol, npvsol(:,:,i-1)./nnsol, npvsol(:,:,i-1));
    A00(1+ndim:2*ndim,1) = reshape(dvv,[],1);
    A00 = A00 + beta*[npusols; npvsols];
    A00 = -A00;

    A2 = zeros(2*ndim,n);
    A2(1:ndim,:) = beta*npsols.*exp(1i * nx_vec * kk) .* (ones(ndim,1) * rk(1,:));
    A2(1+ndim:2*ndim,:) = beta*npsols.*exp(1i * nx_vec * kk) .* (ones(ndim,1) * rk(2,:));
    
    % update the posterior mean and posterior covariance
    mu = mu0 + (a0 + a1 * mu0) * dt + (R0 * A2') * InvBoB * (dnn - (A00 + A2 * mu0) * dt);  % A0 = 0, checked with Nan's book, eqs 8.8 and 8.1, Q.D.
    R = R0 + (a1 * R0 + R0* a1' + Sigma_u*Sigma_u' - (R0*A2') * InvBoB * (R0*A2')')*dt;
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
    plot(dt:dt:N*dt, real(u_post_mean(indd,:)), 'r', 'linewidth',2)
    title(['(a) Gravity mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(indd,:))+2*sqrt(real(u_post_cov(indd,:))), real(u_post_mean(indd,end:-1:1))-2*sqrt(real(u_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = Dim_Ug*2 + 15*(i-1)+4;
    plot(dt:dt:N*dt, real(u_hat(indd,1:N)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean(indd,:)), 'r', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(indd,:))+2*sqrt(real(u_post_cov(indd,:))), real(u_post_mean(indd,end:-1:1))-2*sqrt(real(u_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end

%%
rmsepccPhyDomain

%%
figure
pnx = 32; pny = 32;
[nxx,nyy] = meshgrid(linspace(-pi,pi,pnx), linspace(-pi,pi,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 100;
subplot(2,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('SWE current')

subplot(2,4,5)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('EuDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 1000;
subplot(2,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,6)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('EuDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 10000;
subplot(2,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,7)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('EuDA')
xlabel(['t = ', num2str(dt*ind)])


%--------------------------
ind = 20000;
subplot(2,4,4)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(2,4,8)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(vx.^2 + vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('EuDA')
xlabel(['t = ', num2str(dt*ind)])

%%
%sig = getNoise4EuDA(domain, np, 2*nx, nx, dt, 1)
% sig = getNoise4EuDA(domain, np, 2*nx, nx, dt, 2)
% sig = getNoise4EuDA(domain, np, 2*nx, nx, dt, 3)
% sig = getNoise4EuDA(domain, np, 2*nx, nx, dt, 5)
