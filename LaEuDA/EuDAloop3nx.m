%********************************************************************
%
%  LEMDA: Lagrangian Eulerian Multiscale Data Assimilation (EuDAloop Part)
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

% Eulerian data assimilation with observational data on number density
close all; clc; clear all; 
rng(77); % fix the random number seed to reproduce results
tic 
beta = 1; 
domain = [-pi pi -pi pi];

% generate ocean current
OU_SWEex3 % incompressible flow; only GB modes
timeEuDAocn = toc
save('./uhat/ocn.mat', "u_hat","kk","rk");
% % file_name = sprintf('./data90kkmax4dt4/uhat/ocn.mat');
% % load(file_name) 
nx = 9; ny = nx; ndim = nx^2;
sigma_xy = 0.001; % noise in the Lagrangian tracer equations
sigma_v = 0.00;

npset=(16)*500;
for npi = 1:length(npset)
    np = npset(npi); 

    maxo = solveParticleModelLoop(domain, sigma_xy, sigma_v, np, npi, dt, kk, rk, N, u_hat,beta);
    [npi toc]
end

savedn = 5000; nxset=(1:5)*3+3; np = 8000;
for npi = 1:length(nxset)
    nx = nxset(npi); ny = nx; ndim = nx^2;
    
    hx = (domain(2) - domain(1))/nx;
    hy = (domain(4) - domain(3))/ny;
    npsol = zeros(ny, nx, N);
    npusol = zeros(ny, nx, N);
    npvsol = zeros(ny, nx, N);
     
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
    
    % data processing to fix when one cell does not have floe
    for j=1:N
        [npsol(:,:,j)] = ProcEuNumData(ny, nx, npsol(:,:,j) );
    end
    npusol = npusol*(2*pi)^2/(np*hx*hy);
    npvsol = npvsol*(2*pi)^2/(np*hx*hy);
    npsol = npsol*(2*pi)^2/(np*hx*hy); % number density scaling to around 1.
    timeEuDAdatapro = toc
    
    %save(['/g/data/zv32/seaIceFloeData/cfdata4nx/np' num2str(npi,'%02.f') 'npsol.mat'],"npsol", "npusol", "npvsol");
    
    % smoothing data by moving averages: note smoothing in time doesn't matter;
    % shall not smooth in spacews
    % npusol0 = npusol; npvsol0 = npvsol; npsol0 = npsol;
    % ws = 10;
    % npusol = smoothdata(npusol0,3,"movmean",ws);
    % npvsol = smoothdata(npvsol0,3,"movmean",ws);
    % npsol = smoothdata(npsol0,3,"movmean",ws);
    %npusol = npusol0; npvsol = npvsol0; npsol = npsol0;
    
    
    % mesh
    [nxx,nyy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx));
    nx_vec = [reshape(nxx,[],1), reshape(nyy,[],1)]; % becoming a two column matrix
    
    %% data assimilation
    % sig = getNoise4EuDA(domain, np, 2*nx, nx, dt, 1)
    % sigma_n = sqrt(0.5*sig(1)^2 + 0.5*sig(2)^2) % noise in the number density discretized ODEs
    
    sigma_n =  0.02*npi;
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
    %mu0 = u_hat(:,1); % initial value of posterior mean
    mu0 = 0.001*randn(Dim_U,1)+0.001i*randn(Dim_U,1); 
    n = length(kk(1,:));
    R0 = zeros(n,n); % initial value of posterior covariance
    u_post_mean = zeros(n,N); % posterior mean
    u_post_mean(:,1) = mu0;
    u_post_cov = zeros(n,N); % posterior covariance
    u_post_cov(:,1) = diag(R0); % only save the diagonal elements
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
        
        %     if(mod(i,5000)==0)
        %         figure;
        %         subplot(2,3,1); imagesc(npsol(:,:,i) - npsol(:,:,i-1) ); colorbar; title(['T =', num2str(i*dt)],'fontsize',14)
        %         subplot(2,3,2); imagesc(npusol(:,:,i) - npusol(:,:,i-1)); colorbar
        %         subplot(2,3,3); imagesc(npvsol(:,:,i) - npvsol(:,:,i-1)); colorbar
        %         subplot(2,3,2); imagesc(dvu); colorbar
        %         subplot(2,3,3); imagesc(dvv); colorbar
        %
        %         subplot(2,3,4); imagesc(npsol(:,:,i)); colorbar
        %         subplot(2,3,5); imagesc(npusol(:,:,i)./nnsol ); colorbar
        %         subplot(2,3,6); imagesc(npvsol(:,:,i)./nnsol ); colorbar
        %     end
        
        A2 = zeros(2*ndim,n);
        A2(1:ndim,:) = beta*npsols.*exp(1i * nx_vec * kk) .* (ones(ndim,1) * rk(1,:));
        A2(1+ndim:2*ndim,:) = beta*npsols.*exp(1i * nx_vec * kk) .* (ones(ndim,1) * rk(2,:));
        %A2 = A2;
        
        %     if(i>10000)
        %         tem = 1 + 1*(i-10000)*dt;
        %         A2 = A2/tem;
        %         dnn = dnn/tem;
        %         A00 = A00/tem;
        %     end
        
        % update the posterior mean and posterior covariance
        mu = mu0 + (a0 + a1 * mu0) * dt + (R0 * A2') * InvBoB * (dnn - (A00 + A2 * mu0) * dt);  % A0 = 0, checked with Nan's book, eqs 8.8 and 8.1, Q.D.
        R = R0 + (a1 * R0 + R0* a1' + Sigma_u*Sigma_u' - (R0*A2') * InvBoB * (R0*A2')')*dt;
        u_post_mean(:,i) = mu;
        u_post_cov(:,i) = diag(R);
        mu0 = mu;
        R0 = R;
        %     if (mod(i,10000)==0)
        %         mu0 = u_hat(:,i);
        %         R0 = 0.0*R;
        %     end
    end
    timeEuDA = toc
    %save(['/g/data/zv32/seaIceFloeData/cfdata4nx/np' num2str(npi,'%02.f') 'uhat.mat'],"u_post_mean", "u_post_cov");
    save(['./uhat/EuDAuhat' num2str(npi,'%02.f') 'nxx.mat'],"u_post_mean", "u_post_cov");


    
    
    %% The following lines are for plotting the results
    %u_post_mean = abs(u_post_mean);
    figure
    for i = 1:4
        subplot(4,2,2*i-1)
        hold on
        %indd = mod(24*(i-1)+1,40); % for kmax = 3
        indd = Dim_Ug*2 + 13*(i-1)+1;
        plot(dt:dt:N*dt, real(u_hat(indd,1:N)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean(indd,:)), 'r', 'linewidth',2)
        title(['(a) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
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
    ind = 1000;
    subplot(2,4,1)
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
    set(gca,'fontsize',16)
    box on
    ylabel('SWE current')
    
    subplot(2,4,5)
    vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
    vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
    vx = reshape(real(vx), pny, pnx);
    vy = reshape(real(vy), pny, pnx);
    vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
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
    subplot(2,4,2)
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
    set(gca,'fontsize',16)
    box on
    %ylabel('SWE current')
    
    subplot(2,4,6)
    vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
    vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
    vx = reshape(real(vx), pny, pnx);
    vy = reshape(real(vy), pny, pnx);
    vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
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
    ind = 30000;
    subplot(2,4,3)
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
    set(gca,'fontsize',16)
    box on
    %ylabel('SWE current')
    
    subplot(2,4,7)
    vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
    vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
    vx = reshape(real(vx), pny, pnx);
    vy = reshape(real(vy), pny, pnx);
    vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
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
    ind = 50000;
    subplot(2,4,4)
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
    set(gca,'fontsize',16)
    box on
    %ylabel('SWE current')
    
    subplot(2,4,8)
    vx = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(1,:)));
    vy = exp(1i * nx_vec * kk) * (u_post_mean(:,ind) .* transpose(rk(2,:)));
    vx = reshape(real(vx), pny, pnx);
    vy = reshape(real(vy), pny, pnx);
    vx = real(vx); vy = real(vy); vc = sqrt(0.5*vx.^2 + 0.5*vy.^2);
    hold on; contourf(nxx,nyy,vc,40,'edgecolor','none')
    colorbar
    quiver(xx, yy, vx, vy, 'linewidth',1.5)
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    set(gca,'fontsize',16)
    box on
    ylabel('EuDA')
    xlabel(['t = ', num2str(dt*ind)])
    
end

%%
errpcc = zeros(2, length(nxset));
for npi = 1:length(nxset)
    
    load(['./err/err' num2str(npi,'%01.f') 'nxx.mat'], "err")
    
    errpcc(1, npi) = mean(0.5*(err(20000:end, 1) + err(20000:end, 2) ) );
    errpcc(2, npi) = mean(0.5*(err(20000:end, 3) + err(20000:end, 4) ) );
end

%%
figure
subplot(2,3,3)
hold on
plot(nxset, errpcc(1,:), '-*b', 'linewidth',2)
title('Normalised RMSE','fontsize',14)
set(gca,'fontsize',16)
%legend('vx','vy')
box on
xlabel('Grid size')

subplot(2,3,6)
hold on
plot(nxset, errpcc(2,:), '-*b', 'linewidth',2)
title('PCC','fontsize',14)
set(gca,'fontsize',16)
%legend('vx','vy')
box on
xlabel('Grid size')

%save('./ws/EuDAkmax4np8k-1.mat')
