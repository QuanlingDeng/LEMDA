nn = nx*2; ns = nn/nx;
[xx,yy] = meshgrid(linspace(-pi,pi,nn), linspace(-pi,pi,nn));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix

figure
ind = 1000;
subplot(4,4,1)
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

subplot(4,4,5)
u = reshape(kron(nusol(:,:,ind),ones(ns)), nn, nn);
v = reshape(kron(nvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*u.^2 + 0.5*v.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('Mom Eqn Sol')

subplot(4,4,9)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('Ocean')

subplot(4,4,13)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
ylabel('EuDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 5000;
subplot(4,4,2)
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Particle counts')

subplot(4,4,6)
u = reshape(kron(nusol(:,:,ind),ones(ns)), nn, nn);
v = reshape(kron(nvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*u.^2 + 0.5*v.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Mom Eqn Sol')

subplot(4,4,10)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on

subplot(4,4,14)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 9000;
subplot(4,4,3)
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Particle counts')

subplot(4,4,7)
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

subplot(4,4,11)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on

subplot(4,4,15)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])


%--------------------------
ind = 10000;
subplot(4,4,4)
nu = reshape(kron(npusol(:,:,ind),ones(ns)), nn, nn);
nv = reshape(kron(npvsol(:,:,ind),ones(ns)), nn, nn);
hold on; vc = sqrt(0.5*nu.^2 + 0.5*nv.^2); contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
%ylabel('Particle counts')

subplot(4,4,8)
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

subplot(4,4,12)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on

subplot(4,4,16)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])

%%
ind = 1000;
subplot(4,4,13)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])

ind = 5000;
subplot(4,4,14)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])

ind = 9000;
subplot(4,4,15)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])

ind = 10000;
subplot(4,4,16)
vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), nn, nn);
vy = reshape(real(vy), nn, nn);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(xx,yy,vc,40,'edgecolor','none'); colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',16)
box on
xlabel(['t = ', num2str(dt*ind)])

%%
return

LaDA

EuDA

LEMDA


subplot(4,4,5)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])


ind = 100;
subplot(4,4,2)
u = reshape(nusol(:,:,ind), nn, nn);
v = reshape(nvsol(:,:,ind), nn, nn);
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on

subplot(4,4,6)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])


ind = 200;
subplot(4,4,3)
u = reshape(nusol(:,:,ind), nn, nn);
v = reshape(nvsol(:,:,ind), nn, nn);
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on

subplot(4,4,7)
nu = reshape(npusol(:,:,ind), nn, nn);
nv = reshape(npvsol(:,:,ind), nn, nn);
quiver(xx, yy, nu, nv, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on
xlabel(['t = ', num2str(dt*ind)])

ind = 1000;
subplot(4,4,4)
u = reshape(nusol(:,:,ind), nn, nn);
v = reshape(nvsol(:,:,ind), nn, nn);
quiver(xx, yy, u, v, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
box on

subplot(4,4,8)
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

