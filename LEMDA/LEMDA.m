%********************************************************************
%
%  LEMDA: Lagrangian Eulerian Multiscale Data Assimilation
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

close all; clc; clear all; 
% This file documents the algorithms for the LEMDA framework
% A few points are: 
% a) simulation is done in a nondimensionalised and idealised setting to
% demonstrate the LEMDA framework 
% b) data saving and plotting are optional

%% General setting
rng(77); % Fix the random number seed to reproduce results
tic % Start timing
beta = 1; % Ocean ice floe drag coefficient, nondimensionalised
domain = [-pi pi -pi pi]; % Simulation domain

%% Generating the ocean current and ice floe tracer synthetical data
OU_SWE_LEMDA 
save('./uhat/ocn.mat', "u_hat","kk","rk");
timeLaDAocn = toc

sigma_xy = 0.001; % noise in the Lagrangian tracer equations
sigv = 0;
np = 3136; % np is the total number of particles in the; nqq observed
% 
maxo = solveParticleModel(domain, sigma_xy, sigv, np, dt, kk, rk, N, u_hat,beta);
timeLEMDAparticle = toc

%% large scale EuDA part
nx = 7; ny = nx; ndim = nx^2;
hx = (domain(2) - domain(1))/nx;
hy = (domain(4) - domain(3))/ny;
npsol = zeros(ny, nx, N);
npusol = zeros(ny, nx, N);
npvsol = zeros(ny, nx, N);
npi = 1; savedn = 5000; 

for j=1:N/savedn
    file_name = sprintf('./data/np01time%03d.mat', j);
    load(file_name)
    
    for l=1:savedn
        indx = ceil( ( FloeX(:,l) - domain(1) )/hx );
        indy = ceil( ( FloeY(:,l) - domain(3) )/hy );
        
        indt = (j-1)*savedn + l;
        for k=1:np
            npsol(indy(k), indx(k), indt) = npsol(indy(k), indx(k), indt) + 1;
            npusol(indy(k), indx(k), indt) = npusol(indy(k), indx(k), indt) + FloeU(k,l);
            npvsol(indy(k), indx(k), indt) = npvsol(indy(k), indx(k), indt) + FloeV(k,l);
        end
    end
end

% Data processing to fix when one cell does not have floe
for j=1:N
    [npsol(:,:,j)] = ProcEuNumData(ny, nx, npsol(:,:,j) );
end
npusol = npusol*(2*pi)^2/(np*hx*hy);
npvsol = npvsol*(2*pi)^2/(np*hx*hy);
npsol = npsol*(2*pi)^2/(np*hx*hy); % number density scaling to around 1.
save(['./data/LEMDA2EuDAnp' num2str(npi,'%02.f') 'npsol.mat'],"npsol", "npusol", "npvsol");
timeEuDAdatapro = toc


% Generating mesh
[nxx,nyy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx));
nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix

% EuDA using conditional Gaussian
sigma_n =  0.01;
cInvBoB = eye(2*ndim)/sigma_n/sigma_n;
cmu0 = 0.001*randn(cDim_U,1)+0.001i*randn(cDim_U,1); 
cn = length(ckk(1,:));
cR0 = zeros(cn,cn); % initial value of posterior covariance
cu_post_mean = zeros(cn,N); % posterior mean
cu_post_mean(:,1) = cmu0;
cu_post_cov = zeros(cn,N); % posterior covariance
cu_post_cov(:,1) = diag(cR0); % only save the diagonal elements
for i = 2:N
%     sigma_n = 0.005 + 0.03*i*dt;
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
    dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1), npvsol(:,:,i-1), npusol(:,:,i-1)./nnsol);
    %dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1)./nnsol, npvsol(:,:,i-1)./nnsol, npusol(:,:,i-1)); % same as above; tested
    dvu = dvv;
    A00(1:ndim,1) = reshape(dvv,[],1);
    dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1), npvsol(:,:,i-1), npvsol(:,:,i-1)./nnsol);
    %dvv = getDivVV(domain, nx, ny, npusol(:,:,i-1)./nnsol, npvsol(:,:,i-1)./nnsol, npvsol(:,:,i-1));
    A00(1+ndim:2*ndim,1) = reshape(dvv,[],1);
    A00 = A00 + beta*[npusols; npvsols];
    A00 = -A00;

    A2 = zeros(2*ndim,cn);
    A2(1:ndim,:) = beta*npsols.*exp(1i * nx_vec * ckk) .* (ones(ndim,1) * crk(1,:));
    A2(1+ndim:2*ndim,:) = beta*npsols.*exp(1i * nx_vec * ckk) .* (ones(ndim,1) * crk(2,:));

    % update the posterior mean and posterior covariance
    cmu = cmu0 + (ca0 + ca1 * cmu0) * dt + (cR0 * A2') * cInvBoB * (dnn - (A00 + A2 * cmu0) * dt);  % A0 = 0, checked with Nan's book, eqs 8.8 and 8.1, Q.D.
    cR = cR0 + (ca1 * cR0 + cR0* ca1' + cSigma_u*cSigma_u' - (cR0*A2') * cInvBoB * (cR0*A2')')*dt;
    cu_post_mean(:,i) = cmu;
    cu_post_cov(:,i) = diag(cR);
    cmu0 = cmu;
    cR0 = cR;
end


timeEuDA = toc
save(['./uhat/LEMDAcuhat' num2str(npi,'%02.f') '.mat'],"cu_post_mean", "cu_post_cov");
%load(['./uhat/LEMDAcuhat' num2str(npi,'%02.f') '.mat'],"cu_post_mean", "cu_post_cov");

rmsepcceuda4lemda % Evaluating EuDA skill scores (RMSE & PCC)

%% Geerating figures
figure
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    indd = 30+2*i-1; % for kmax = 3
    plot(dt:dt:N*dt, real(cu_post_mean(indd,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_hat(indcmod(indd),1:N)), 'b', 'linewidth',2)
    title(['(a) GB mode ( ', num2str(ckk(1,indd)),' , ', num2str(ckk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(cu_post_mean(indd,:))+2*sqrt(real(cu_post_cov(indd,:))), real(cu_post_mean(indd,end:-1:1))-2*sqrt(real(cu_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = 20+3*i-2; 
    plot(dt:dt:N*dt, real(cu_post_mean(indd,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_hat(indcmod(indd),1:N)), 'b', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(cu_post_mean(indd,:))+2*sqrt(real(cu_post_cov(indd,:))), real(cu_post_mean(indd,end:-1:1))-2*sqrt(real(cu_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(ckk(1,indd)),' , ', num2str(ckk(2,indd)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end
uedatime = toc

%% Lagrangian DA (LaDA) for small scale modes with |k|>4
nx = 7; ny = nx; ndim = nx^2; hx = (domain(2) - domain(1))/nx; hy = (domain(4) - domain(3))/ny;
nqq = 10; savedn = 5000;
xl = zeros(ndim*nqq, N); yl = zeros(ndim*nqq, N);
xr = zeros(ndim*nqq, N+1); yr = zeros(ndim*nqq, N+1);
indl = zeros(ndim*nqq,1); indr = zeros(ndim*nqq,1);
FloeXX = zeros(np, N+1); FloeYY = zeros(np, N+1);

for j=1:N/savedn
    file_name = sprintf('./data/np01time%03d.mat', j);
    load(file_name)
    
    FloeXX(:, (j-1)*savedn+1:j*savedn) = FloeX;
    FloeYY(:, (j-1)*savedn+1:j*savedn) = FloeY;
end

datadn = N;
for j=1:N/datadn
    indeu = zeros(ndim,1);
    for k=1:np
        x0 = FloeXX(k, (j-1)*datadn+1); y0 = FloeYY(k, (j-1)*datadn+1);
        indx = ceil((x0+pi)/hx); indy = ceil((y0+pi)/hy);
        cindx = (indx-0.5)*hx - pi; cindy = (indy-0.5)*hy - pi; % center of the cell (indx, indy)
        dis = sqrt( (x0-cindx)^2 + (y0-cindy)^2);
        ind = (indy - 1)*nx + indx;
        
        if(dis < 0.45*hx && indeu(ind)<nqq)
            indeu(ind) = indeu(ind) + 1;
            indl((ind-1)*nqq+indeu(ind), 1) = k;
        end
    end
    
    if(sum(indeu<nqq)>0)
        disp(j)
        disp(sum(indeu<nqq))
        disp("not found enough observing floes in local cells within 0.25hx")
        return
    end
    
    xl(:, (j-1)*datadn+1:j*datadn) = FloeXX(indl, (j-1)*datadn+1:j*datadn);
    yl(:, (j-1)*datadn+1:j*datadn) = FloeYY(indl, (j-1)*datadn+1:j*datadn);
    
    xr(:, (j-1)*datadn+2:j*datadn+1) = FloeXX(indl, (j-1)*datadn+2:j*datadn+1);
    yr(:, (j-1)*datadn+2:j*datadn+1) = FloeYY(indl, (j-1)*datadn+2:j*datadn+1);
end
toc


% Run LaDA in each cell using conditional Gaussian
sigma_v = 0.02;
fn = length(fkk(1,:)) + 2*nqq; %
fInvBoB = eye(2*nqq)/sigma_xy/sigma_xy;
sfu_post_mean = zeros(ndim,fn,N); % total of post_mean
sfu_post_cov = zeros(ndim, fn,N); 
for jy = 1:ny
    jy
    toc
    for jx = 1:nx
        ind = (jy - 1)*nx + jx; indl = (ind-1)*nqq+1; indr = ind*nqq;
        xxl = xl(indl:indr,:); xxr = xr(indl:indr,:); yyl = yl(indl:indr,:); yyr = yr(indl:indr,:);
        
        fmu0 = 0.00*randn(fn,1)+0.0*0.1i*randn(fn,1);
        fmu0(1:2*nqq) = 0;
        fR0 = zeros(fn,fn); % initial value of posterior covariance
        fu_post_mean = zeros(fn,N); % posterior mean
        fu_post_mean(:,1) = fmu0;
        fu_post_cov = zeros(fn,N); % posterior covariance
        fu_post_cov(:,1) = diag(fR0); % only save the diagonal elements
        fA1 = [eye(2*nqq), zeros(2*nqq, fDim_U)];
        fQ = zeros(2*nqq,fDim_U); fQQ = zeros(fDim_U, 2*nqq);
        cU = zeros(2*nqq,1);
        fsa0 = [zeros(2*nqq,1); fa0]; % sum a0 part
        fb2b = fSigma_u*fSigma_u';
        fb2dotb2 = [eye(2*nqq)*sigma_v*sigma_v fQ; fQQ fb2b];

        for i = 2:N
            % matrix for filtering
            xdiff = [xxr(:,i)-xxl(:,i-1); yyr(:,i)-yyl(:,i-1)];
            xdiff(xdiff>pi) = xdiff(xdiff>pi) - 2*pi; % periodic boundary condition
            xdiff(xdiff<-pi) = xdiff(xdiff<-pi) + 2*pi;
            
            fQ(1:nqq,:)      = exp(1i * xxl(:,i-1) * fkk(1,:) + 1i * yyl(:,i-1) * fkk(2,:)) .* (ones(nqq,1) * frk(1,:));
            fQ(nqq+1:2*nqq,:) = exp(1i * xxl(:,i-1) * fkk(1,:) + 1i * yyl(:,i-1) * fkk(2,:)) .* (ones(nqq,1) * frk(2,:));
            fQ = beta*fQ;
            fsa1 = [-beta*eye(2*nqq) fQ; fQQ fa1];
            
            nx_vec = [xxl(:,i-1) yyl(:,i-1)]; 
            cU(1:nqq) = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(1,:)));
            cU(nqq+1:2*nqq) = exp(1i * nx_vec * ckk) * (cu_post_mean(:,ind) .* transpose(crk(2,:)));
            cU = beta*cU;
            fsa0 = [cU; fa0];
            
            % update the posterior mean and posterior covariance
            fmu = fmu0 + (fsa0 + fsa1 * fmu0) * dt + (fR0 * fA1') * fInvBoB * (xdiff - fA1 * fmu0 * dt);  % A0 = 0
            fR = fR0 + (fsa1 * fR0 + fR0* fsa1' + fb2dotb2 - (fR0*fA1') * fInvBoB * (fR0*fA1')')*dt;
            fu_post_mean(:,i) = fmu;
            fu_post_cov(:,i) = diag(real(fR));
            
            fmu0 = fmu;
            fR0 = fR;
        end
        
        sfu_post_mean(ind, :, :) = fu_post_mean;
        sfu_post_cov(ind, :, :) = fu_post_cov;
    end
end
lamdatime = toc
save(['./uhat/LEMDAfuhat' num2str(npi,'%02.f') '.mat'],'-v7.3',"sfu_post_mean", "sfu_post_cov");

%% evaluating the LEMDA skill scores
rmsepcc4lemda

%%
figure
indcell = 11;
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    indd = 2*nqq+5+2*i; %2*nqq+mod(2*(i-1)+1,79); % for kmax = 3
    plot(dt:dt:N*dt, real(u_hat(indfmod(indd-2*nqq),1:N)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, reshape(real(sfu_post_mean(indcell,indd,:)),[],1), 'r', 'linewidth',2)
    title(['(a) GB mode ( ', num2str(fkk(1,indd-2*nqq)),' , ', num2str(fkk(2,indd-2*nqq)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [reshape(real(sfu_post_mean(indcell,indd,:))+2*sqrt(real(sfu_post_cov(indcell,indd,:))),[],1); reshape(real(sfu_post_mean(indcell,indd,end:-1:1))-2*sqrt(real(sfu_post_cov(indcell,indd,end:-1:1))),[],1)]','r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = 2*nqq+20+2*i; %2*nqq+mod(3*(i-1)+1,79); 
    plot(dt:dt:N*dt, real(u_hat(indfmod(indd-2*nqq),1:N)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, reshape(real(sfu_post_mean(indcell,indd,:)),[],1), 'r', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [reshape(real(sfu_post_mean(indcell,indd,:))+2*sqrt(real(sfu_post_cov(indcell,indd,:))),[],1); reshape(real(sfu_post_mean(indcell,indd,end:-1:1))-2*sqrt(real(sfu_post_cov(indcell,indd,end:-1:1))),[],1)]','r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(fkk(1,indd-2*nqq)),' , ', num2str(fkk(2,indd-2*nqq)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end


%% Plotting figures
figure
pnx = 35; pny = 35; h = 2*pi/pnx; fny = pny/ny; fnx = fny;
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
fvx = 0*vx; fvy = fvx;
for jy = 1:ny
    for jx = 1:nx
        indu = (jy - 1)*nx + jx;
        utem = reshape(sfu_post_mean(indu, 2*nqq+1:end, ind),[],1);
        
        temfvx = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(1,:)));
        temfvy = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(2,:)));
        temfvx = reshape(real(temfvx), pny, pnx);
        temfvy = reshape(real(temfvy), pny, pnx);
        
        fvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
        fvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
    end
end
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
fvx = 0*vx; fvy = fvx;
for jy = 1:ny
    for jx = 1:nx
        indu = (jy - 1)*nx + jx;
        utem = reshape(sfu_post_mean(indu, 2*nqq+1:end, ind),[],1);
        
        temfvx = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(1,:)));
        temfvy = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(2,:)));
        temfvx = reshape(real(temfvx), pny, pnx);
        temfvy = reshape(real(temfvy), pny, pnx);
        
        fvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
        fvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
    end
end
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
fvx = 0*vx; fvy = fvx;
for jy = 1:ny
    for jx = 1:nx
        indu = (jy - 1)*nx + jx;
        utem = reshape(sfu_post_mean(indu, 2*nqq+1:end, ind),[],1);
        
        temfvx = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(1,:)));
        temfvy = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(2,:)));
        temfvx = reshape(real(temfvx), pny, pnx);
        temfvy = reshape(real(temfvy), pny, pnx);
        
        fvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
        fvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
    end
end
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


% save('./ws/LEMDA2.mat','-v7.3')