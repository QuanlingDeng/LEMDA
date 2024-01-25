%% Lagragian DA part
%nqq = nqq*ndim;
sigma_v = 0.02;
fn = length(fkk(1,:)) + 2*nqq; %
fInvBoB = eye(2*nqq)/sigma_xy/sigma_xy;
fmu0 = 0.00*randn(fn,1)+0.0*0.1i*randn(fn,1);
fmu0(1:2*nqq) = 0;
fR0 = zeros(fn,fn); % initial value of posterior covariance
ffu_post_mean = zeros(fn,N); % posterior mean
ffu_post_mean(:,1) = fmu0;
ffu_post_cov = zeros(fn,N); % posterior covariance
ffu_post_cov(:,1) = diag(fR0); % only save the diagonal elements
fA1 = [eye(2*nqq), zeros(2*nqq, fDim_U)];
fQ = zeros(2*nqq,fDim_U); fQQ = zeros(fDim_U, 2*nqq);
cU = zeros(2*nqq,1);
fsa0 = [zeros(2*nqq,1); fa0]; % sum a0 part
fb2b = fSigma_u*fSigma_u';
fb2dotb2 = [eye(2*nqq)*sigma_v*sigma_v fQ; fQQ fb2b];

for i = 2:N
    % matrix for filtering
    xdiff = [xr(:,i)-xl(:,i-1); yr(:,i)-yl(:,i-1)];
    xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
    xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;
    
    fQ(1:nqq,:)      = exp(1i * xl(:,i-1) * fkk(1,:) + 1i * yl(:,i-1) * fkk(2,:)) .* (ones(nqq,1) * frk(1,:));
    fQ(nqq+1:2*nqq,:) = exp(1i * xl(:,i-1) * fkk(1,:) + 1i * yl(:,i-1) * fkk(2,:)) .* (ones(nqq,1) * frk(2,:));
    fQ = beta*fQ;
    fsa1 = [-beta*eye(2*nqq) fQ; fQQ fa1];
    
    nx_vec = [xl(:,i-1) yl(:,i-1)];
    cU(1:nqq) = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
    cU(nqq+1:2*nqq) = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
    cU = beta*cU;
    fsa0 = [cU; fa0];
    
    % update the posterior mean and posterior covariance
    fmu = fmu0 + (fsa0 + fsa1 * fmu0) * dt + (fR0 * fA1') * fInvBoB * (xdiff - fA1 * fmu0 * dt);  % A0 = 0
    fR = fR0 + (fsa1 * fR0 + fR0* fsa1' + fb2dotb2 - (fR0*fA1') * fInvBoB * (fR0*fA1')')*dt;
    ffu_post_mean(:,i) = fmu;
    ffu_post_cov(:,i) = diag(real(fR));
    
    fmu0 = fmu;
    fR0 = fR;
end

%
figure
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    indd = 2*nqq+28+2*i-1; % for kmax = 3
    plot(dt:dt:N*dt, real(ffu_post_mean(indd,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_hat(indfmod(indd-2*nqq),1:N)), 'b', 'linewidth',2)
    title(['(a) GB mode ( ', num2str(fkk(1,indd-2*nqq)),' , ', num2str(fkk(2,indd-2*nqq)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(ffu_post_mean(indd,:))+2*sqrt(real(ffu_post_cov(indd,:))), real(ffu_post_mean(indd,end:-1:1))-2*sqrt(real(ffu_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = 2*nqq+19+2*i-1; 
    plot(dt:dt:N*dt, real(ffu_post_mean(indd,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_hat(indfmod(indd-2*nqq),1:N)), 'b', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(ffu_post_mean(indd,:))+2*sqrt(real(ffu_post_cov(indd,:))), real(ffu_post_mean(indd,end:-1:1))-2*sqrt(real(ffu_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(fkk(1,indd-2*nqq)),' , ', num2str(fkk(2,indd-2*nqq)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end


% plots
figure
pnx = 28; pny = 28; h = 2*pi/pnx; fny = pny/ny; fnx = fny;
[nxx,nyy] = meshgrid(linspace(-pi+0.5*h,pi-0.5*h,pnx), linspace(-pi+0.5*h,pi-0.5*h,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 10000;
subplot(4,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('Ocean current')

subplot(4,4,5)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('EuDA part')
%xlabel(['t = ', num2str(dt*ind)])

vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
temfvx = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(frk(1,:)));
temfvy = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(frk(2,:)));
fvx = reshape(real(temfvx), pny, pnx);
fvy = reshape(real(temfvy), pny, pnx);
subplot(4,4,9)
fvx = real(fvx); fvy = real(fvy); fvc = sqrt(0.5*fvx.^2 + 0.5*fvy.^2);
hold on; contourf(nxx,nyy,fvc,40,'edgecolor','none')
colorbar
quiver(xx, yy, fvx, fvy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('LEMDA part')

subplot(4,4,13)
vx = vx + fvx; vy = vy + fvy;    
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 30000;
subplot(4,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('Ocean current')

subplot(4,4,6)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('EuDA')
%xlabel(['t = ', num2str(dt*ind)])

vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
temfvx = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(frk(1,:)));
temfvy = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(frk(2,:)));
fvx = reshape(real(temfvx), pny, pnx);
fvy = reshape(real(temfvy), pny, pnx);
subplot(4,4,10)
fvx = real(fvx); fvy = real(fvy); fvc = sqrt(0.5*fvx.^2 + 0.5*fvy.^2);
hold on; contourf(nxx,nyy,fvc,40,'edgecolor','none')
colorbar
quiver(xx, yy, fvx, fvy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on

subplot(4,4,14)
vx = vx + fvx; vy = vy + fvy;    
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 50000;
subplot(4,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('Ocean current')

subplot(4,4,7)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('EuDA')
%xlabel(['t = ', num2str(dt*ind)])

vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
temfvx = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(frk(1,:)));
temfvy = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(frk(2,:)));
fvx = reshape(real(temfvx), pny, pnx);
fvy = reshape(real(temfvy), pny, pnx);

subplot(4,4,11)
fvx = real(fvx); fvy = real(fvy); fvc = sqrt(0.5*fvx.^2 + 0.5*fvy.^2);
hold on; contourf(nxx,nyy,fvc,40,'edgecolor','none')
colorbar
quiver(xx, yy, fvx, fvy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on

subplot(4,4,15)
vx = vx + fvx; vy = vy + fvy;    
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
toc
xlabel(['t = ', num2str(dt*ind)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% kk %%%%%%%%%%%%%%%%


%% Lagragian DA part
%nqq = nqq*ndim;
fn = length(kk(1,:)) + 2*nqq; %
fInvBoB = eye(2*nqq)/sigma_xy/sigma_xy;
fmu0 = 0.00*randn(fn,1)+0.0*0.1i*randn(fn,1);
fmu0(1:2*nqq) = 0;
fR0 = zeros(fn,fn); % initial value of posterior covariance
ffu_post_mean = zeros(fn,N); % posterior mean
ffu_post_mean(:,1) = fmu0;
ffu_post_cov = zeros(fn,N); % posterior covariance
ffu_post_cov(:,1) = diag(fR0); % only save the diagonal elements
fA1 = [eye(2*nqq), zeros(2*nqq, Dim_U)];
fQ = zeros(2*nqq,Dim_U); fQQ = zeros(Dim_U, 2*nqq);
cU = zeros(2*nqq,1);
fsa0 = [zeros(2*nqq,1); a0]; % sum a0 part
fb2b = Sigma_u*Sigma_u';
fb2dotb2 = [eye(2*nqq)*sigma_v*sigma_v fQ; fQQ fb2b];

for i = 2:N
    % matrix for filtering
    xdiff = [xr(:,i)-xl(:,i-1); yr(:,i)-yl(:,i-1)];
    xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
    xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;
    
    fQ(1:nqq,:)      = exp(1i * xl(:,i-1) * kk(1,:) + 1i * yl(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(1,:));
    fQ(nqq+1:2*nqq,:) = exp(1i * xl(:,i-1) * kk(1,:) + 1i * yl(:,i-1) * kk(2,:)) .* (ones(nqq,1) * rk(2,:));
    fQ = beta*fQ;
    fsa1 = [-beta*eye(2*nqq) fQ; fQQ a1];
    
    nx_vec = [xl(:,i-1) yl(:,i-1)];
    cU(1:nqq) = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
    cU(nqq+1:2*nqq) = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
    cU = beta*cU;
    fsa0 = [0*cU; a0];
    
    % update the posterior mean and posterior covariance
    fmu = fmu0 + (fsa0 + fsa1 * fmu0) * dt + (fR0 * fA1') * fInvBoB * (xdiff - fA1 * fmu0 * dt);  % A0 = 0
    fR = fR0 + (fsa1 * fR0 + fR0* fsa1' + fb2dotb2 - (fR0*fA1') * fInvBoB * (fR0*fA1')')*dt;
    ffu_post_mean(:,i) = fmu;
    ffu_post_cov(:,i) = diag(real(fR));
    
    fmu0 = fmu;
    fR0 = fR;
end


%% plots
figure
pnx = 28; pny = 28; h = 2*pi/pnx; fny = pny/ny; fnx = fny;
[nxx,nyy] = meshgrid(linspace(-pi+0.5*h,pi-0.5*h,pnx), linspace(-pi+0.5*h,pi-0.5*h,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 10000;
subplot(4,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('Ocean current')

subplot(4,4,5)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('EuDA part')
%xlabel(['t = ', num2str(dt*ind)])

vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
temfvx = exp(1i * nx_vec * kk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
temfvy = exp(1i * nx_vec * kk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
fvx = reshape(real(temfvx), pny, pnx);
fvy = reshape(real(temfvy), pny, pnx);
subplot(4,4,9)
fvx = real(fvx); fvy = real(fvy); fvc = sqrt(0.5*fvx.^2 + 0.5*fvy.^2);
hold on; contourf(nxx,nyy,fvc,40,'edgecolor','none')
colorbar
quiver(xx, yy, fvx, fvy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('LEMDA part')

subplot(4,4,13)
vx = vx + fvx; vy = vy + fvy;    
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 30000;
subplot(4,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('Ocean current')

subplot(4,4,6)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('EuDA')
%xlabel(['t = ', num2str(dt*ind)])

vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
temfvx = exp(1i * nx_vec * kk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
temfvy = exp(1i * nx_vec * kk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
fvx = reshape(real(temfvx), pny, pnx);
fvy = reshape(real(temfvy), pny, pnx);
subplot(4,4,10)
fvx = real(fvx); fvy = real(fvy); fvc = sqrt(0.5*fvx.^2 + 0.5*fvy.^2);
hold on; contourf(nxx,nyy,fvc,40,'edgecolor','none')
colorbar
quiver(xx, yy, fvx, fvy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on

subplot(4,4,14)
vx = vx + fvx; vy = vy + fvy;    
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('LEMDA')
xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 50000;
subplot(4,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('Ocean current')

subplot(4,4,7)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('EuDA')
%xlabel(['t = ', num2str(dt*ind)])

vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
temfvx = exp(1i * nx_vec * kk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(rk(1,:)));
temfvy = exp(1i * nx_vec * kk) * (ffu_post_mean(2*nqq+1:end,ind) .* transpose(rk(2,:)));
fvx = reshape(real(temfvx), pny, pnx);
fvy = reshape(real(temfvy), pny, pnx);

subplot(4,4,11)
fvx = real(fvx); fvy = real(fvy); fvc = sqrt(0.5*fvx.^2 + 0.5*fvy.^2);
hold on; contourf(nxx,nyy,fvc,40,'edgecolor','none')
colorbar
quiver(xx, yy, fvx, fvy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on

subplot(4,4,15)
vx = vx + fvx; vy = vy + fvy;    
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
toc
xlabel(['t = ', num2str(dt*ind)])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure
pnx = 32; pny = 32;
[nxx,nyy] = meshgrid(linspace(-pi,pi,pnx), linspace(-pi,pi,pny));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

xx = nxx; yy = nyy;
ind = 10000;
subplot(4,4,1)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('Ocean current')

subplot(4,4,9)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('EuDA')
%xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 30000;
subplot(4,4,2)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('Ocean current')

subplot(4,4,10)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('EuDA')
%xlabel(['t = ', num2str(dt*ind)])

%--------------------------
ind = 50000;
subplot(4,4,3)
vx = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(1,:)));
vy = exp(1i * nx_vec * kk) * (u_hat(:,ind) .* transpose(rk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
%ylabel('Ocean current')

subplot(4,4,11)
vx = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
vy = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
vx = reshape(real(vx), pny, pnx);
vy = reshape(real(vy), pny, pnx);
vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
colorbar
quiver(xx, yy, vx, vy, 'linewidth',1.5)
xlim([-pi, pi ])
ylim([-pi, pi ])
set(gca,'fontsize',24); set(gca,'linewidth',2)
box on
ylabel('EuDA')
xlabel(['t = ', num2str(dt*ind)])