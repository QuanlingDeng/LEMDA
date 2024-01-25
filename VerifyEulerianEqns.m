% Verify the observational data on number density by comparing CL solver
% and particles
rng(77); % fix the random number seed to reproduce results
tic
beta = 0.1; 
domain = [-pi pi -pi pi];

% generate ocean current
%OU_SWEex2 % incompressible flow; only GB modes
nt = N; toc

nx = 15; 1*(2 * K_max + 1); ny = nx; ndim = nx^2;
% 
% %% get the number density data; 230400=480^2; 129600=360^2
% sigma_xy = 0.001; % noise in the Lagrangian tracer equations
% sigma_v = 0.01;
np = 2916; % nqq = 80; % np is the total number of particles in the; nqq observed
% % 
% maxo = solveParticleModel(domain, sigma_xy, sigma_v, np, dt, kk, rk, N, u_hat,beta);
% toc


%%
hx = (domain(2) - domain(1))/nx;
hy = (domain(4) - domain(3))/ny;
nsol = zeros(ny, nx, nt);
npsol = zeros(ny, nx, nt);
nusol = zeros(ny, nx, nt);
nvsol = zeros(ny, nx, nt);
npusol = zeros(ny, nx, nt);
npvsol = zeros(ny, nx, nt);
for j=1:nt
    file_name = sprintf('./data/time%05d.mat', j);
    load(file_name)   
    
    indx = ceil( ( La(:,1) - domain(1) )/hx );
    indy = ceil( ( La(:,2)  - domain(3) )/hy );

    for k=1:np
        npsol(indy(k), indx(k), j) = npsol(indy(k), indx(k), j) + 1;
        npusol(indy(k), indx(k), j) = npusol(indy(k), indx(k), j) + La(k,3);
        npvsol(indy(k), indx(k), j) = npvsol(indy(k), indx(k), j) + La(k,4);
    end
end
npusol = npusol*(2*pi)^2/(np*hx*hy);
npvsol = npvsol*(2*pi)^2/(np*hx*hy);
npsol = npsol*(2*pi)^2/(np*hx*hy);
toc

% data processing to fix when one cell does not have floe
for j=1:nt
    [npsol(:,:,j)] = ProcEuNumData(ny, nx, npsol(:,:,j) );
end

sol = npsol(:,:,1);
usol = npusol(:,:,1); %.*rho(:,:,1);
vsol = npvsol(:,:,1); %.*rho(:,:,1);
nsol(:,:,1) = sol;
nusol(:,:,1) = real(usol);
nvsol(:,:,1) = real(vsol);
[xx,yy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix
for i = 2:nt % generating the number density data
    uox = exp(1i * nx_vec * kk) * (u_hat(:,i-1) .* transpose(rk(1,:)));
    uoy = exp(1i * nx_vec * kk) * (u_hat(:,i-1) .* transpose(rk(2,:)));
    uox = reshape(real(uox), ny, nx);
    uoy = reshape(real(uoy), ny, nx);
    vx = usol./sol; vy = vsol./sol;
    
    forcex = beta*(uox.*npsol(:,:,i) - npusol(:,:,i));
    forcey = beta*(uoy.*npsol(:,:,i) - npvsol(:,:,i));

    ssol = solveNumDensity(domain, nx, ny, dt, vx, vy, sol*0.0, sol, 0.0);
    nsol(:,:,i) = ssol; sol = ssol;
    
    ussol = solveNumDensity(domain, nx, ny, dt, vx, vy, forcex, usol, 0.0);
    nusol(:,:,i) = ussol; usol = ussol;
    
    vssol = solveNumDensity(domain, nx, ny, dt, vx, vy, forcey, vsol, 0.0);
    nvsol(:,:,i) = vssol; vsol = vssol;
end

%%
% npusol = npusol./npsol; npvsol = npvsol./npsol;
% nusol = nusol./nsol; nvsol = nvsol./nsol;

nn = nx*2; ns = nn/nx;
[xx,yy] = meshgrid(linspace(-pi,pi,nn), linspace(-pi,pi,nn));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix

figure
ind = 1000;
subplot(3,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('SWE current')

subplot(3,4,5)
% nu = reshape(npusol(:,:,ind), nn, nn);
% nv = reshape(npvsol(:,:,ind), nn, nn);
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('Particle counts')

subplot(3,4,9)
u = reshape(kron(nusol(:,:,ind),ones(ns)), nn, nn);
v = reshape(kron(nvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*u.^2 + 0.5*v.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('Mom Eqn Sol')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 5000;
subplot(3,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(3,4,6)
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Particle counts')

subplot(3,4,10)
u = reshape(kron(nusol(:,:,ind),ones(ns)), nn, nn);
v = reshape(kron(nvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*u.^2 + 0.5*v.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Mom Eqn Sol')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 10000;
subplot(3,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(3,4,7)
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Particle counts')

subplot(3,4,11)
u = reshape(kron(nusol(:,:,ind),ones(ns)), nn, nn);
v = reshape(kron(nvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*u.^2 + 0.5*v.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Mom Eqn Sol')
xlabel(['t = ', num2str(dt*ind)])


%--------------------------
ind = 39900;
subplot(3,4,4)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind+1) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
hold on; vc = sqrt(0.5*vx.^2 + 0.5*vy.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('SWE current')

subplot(3,4,8)
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Particle counts')

subplot(3,4,12)
u = reshape(kron(nusol(:,:,ind),ones(ns)), nn, nn);
v = reshape(kron(nvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*u.^2 + 0.5*v.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Mom Eqn Sol')
xlabel(['t = ', num2str(dt*ind)])

%%
return

LaDA

EuDA

LEMDA


subplot(3,4,5)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])


ind = 100;
subplot(3,4,2)
u = reshape(nusol(:,:,ind), nn, nn);
v = reshape(nvsol(:,:,ind), nn, nn);
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on

subplot(3,4,6)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])


ind = 200;
subplot(3,4,3)
u = reshape(nusol(:,:,ind), nn, nn);
v = reshape(nvsol(:,:,ind), nn, nn);
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on

subplot(3,4,7)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])

ind = 1000;
subplot(3,4,4)
u = reshape(nusol(:,:,ind), nn, nn);
v = reshape(nvsol(:,:,ind), nn, nn);
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on

subplot(3,4,8)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])


return
%% number density
figure
nsol = real(nsol);
ind = 1;
subplot(2,4,1)
imagesc(nsol(:,:,ind))
ylabel('ntum ntumberDensity')
set(gca,'fontsize',16)
xlabel(['t = ', num2str(dt*ind)])
subplot(2,4,5)
imagesc(npsol(:,:,ind))
ylabel('Particle counts')
set(gca,'fontsize',16)

ind = 10;
subplot(2,4,2)
imagesc(nsol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
set(gca,'fontsize',16)
subplot(2,4,6)
imagesc(npsol(:,:,ind))
set(gca,'fontsize',16)

ind = 50;
subplot(2,4,3)
imagesc(nsol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
set(gca,'fontsize',16)
subplot(2,4,7)
imagesc(npsol(:,:,ind))
set(gca,'fontsize',16)

ind = 100;
subplot(2,4,4)
imagesc(nsol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
set(gca,'fontsize',16)
subplot(2,4,8)
imagesc(npsol(:,:,ind))
set(gca,'fontsize',16)

%% momemtum in x
figure
nusol = real(nusol);
ind = 1;
subplot(2,4,1)
imagesc(nusol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
subplot(2,4,5)
imagesc(npusol(:,:,ind))


ind = 1000;
subplot(2,4,2)
imagesc(nusol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
subplot(2,4,6)
imagesc(npusol(:,:,ind))

ind = 2000;
subplot(2,4,3)
imagesc(nusol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
subplot(2,4,7)
imagesc(npusol(:,:,ind))

ind = 3000;
subplot(2,4,4)
imagesc(nusol(:,:,ind))
xlabel(['t = ', num2str(dt*ind)])
subplot(2,4,8)
imagesc(npusol(:,:,ind))

