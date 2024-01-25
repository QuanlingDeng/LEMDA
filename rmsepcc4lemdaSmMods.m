%% RMSE root mean square error and PCC pattern cross correlation in physical domain
rmse = zeros(N,2); rrmse = zeros(N,2); pcc = zeros(N,2); 

% [xx,yy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx)); % 2nx or nx; no difference
% nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix
pnx = 3*ny; pny = 3*nx; fny = 3; fnx = 3; h = 2*pi/pny;
[xx,yy] = meshgrid(linspace(-pi+0.5*h,pi-0.5*h,pnx), linspace(-pi+0.5*h,pi-0.5*h,pny));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix

ffQ = zeros(2*pnx*pny, fDim_U); fuhat = zeros(fDim_U, N);

for j=2:N
    vx = real(u_hat(indfmod,j)); vy = imag(u_hat(indfmod,j));  
    
    fvx = zeros(pny); fvy = zeros(pny);
    for jy = 1:ny
        for jx = 1:nx
            indu = (jy - 1)*nx + jx;
            utem = reshape(sfu_post_mean(indu, 2*nqq+1:end, j),[],1);
            
            temfvx = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(1,:)));
            temfvy = exp(1i * nx_vec * fkk) * (utem .* transpose(frk(2,:)));
            temfvx = reshape(real(temfvx), pny, pnx);
            temfvy = reshape(real(temfvy), pny, pnx);

            fvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvx((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
            fvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx) = temfvy((jy-1)*fny+1:jy*fny, (jx-1)*fnx+1:jx*fnx);
        end
    end
    fvx = reshape(fvx,[],1); fvy = reshape(fvy,[],1);
    
    ffQ(1:end/2,:)     = exp(1i * nx_vec * fkk) .* (ones(pnx*pny,1) * frk(1,:));
    ffQ(end/2+1:end,:) = exp(1i * nx_vec * fkk) .* (ones(pnx*pny,1) * frk(2,:));
    fvv = [fvx; fvy];
    %vvda = ffQ\fvv; 
    vvda = lscov(ffQ, fvv); fuhat(:,j) = vvda;
    vxda = real(vvda); vyda = imag(vvda); 
       
    rmse(j,1) = mean( (vx - vxda).^2 ); 
    rrmse(j,1) = mean( vx.^2 ); 
    mvx = mean(vx); mvxda = mean(vxda);
    pcc(j,1) = dot(vx - mvx, vxda -mvxda ) / sqrt( dot(vx-mvx,vx-mvx) * dot(vxda-mvxda,vxda-mvxda) );

    rmse(j,2) = mean( (vy - vyda).^2 ); 
    rrmse(j,2) = mean( vy.^2 ); 
    mvy = mean(vy); mvyda = mean(vyda);
    pcc(j,2) = dot(vy-mvy, vyda-mvyda) / sqrt( dot(vy-mvy,vy-mvy) * dot(vyda-mvyda,vyda-mvyda) );
end

rrmse = sqrt(rmse./rrmse); %sqrt(rmse./mean(rrmse(end/2:end)) );
rmse = sqrt(rmse); 

err = [rrmse pcc];
%save(['./err/nsv' num2str(nsv,'%01.f') 'err.mat'], "err")
save(['./uhat/LEMDAnpi' num2str(npi,'%01.f') 'SmallMods.mat'], "fuhat")
save(['./err/LEMDAnpi' num2str(npi,'%01.f') 'SmallModsErr.mat'], "err")

%%
figure
subplot(2,3,1)
hold on
% plot(dt:dt:N*dt, rrmse(:,1), ':', 'linewidth',2)
% plot(dt:dt:N*dt, rrmse(:,2), '--', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*rrmse(:,1).^2+0.5*rrmse(:,2).^2), '-', 'linewidth',2)
title('Normalised RMSE','fontsize',24)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on
xlabel('t')

subplot(2,3,4)
hold on
% plot(dt:dt:N*dt, pcc(:,1), ':', 'linewidth',2)
% plot(dt:dt:N*dt, pcc(:,2), '--', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*pcc(:,1).^2+0.5*pcc(:,2).^2), '-', 'linewidth',2)
title('PCC','fontsize',24)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on
xlabel('t')


%%

subplot(4,4,16)
hold on
% plot(dt:dt:N*dt, pcc(:,1), ':', 'linewidth',2)
% plot(dt:dt:N*dt, pcc(:,2), '--', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*rrmse(:,1).^2+0.5*rrmse(:,2).^2), '-', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*pcc(:,1).^2+0.5*pcc(:,2).^2), '-', 'linewidth',2)
title('PCC & RMSE','fontsize',24)
set(gca,'fontsize',24); set(gca,'linewidth',2)
%legend('vx','vy')
box on
xlabel('t')

%%
figure
indcell = 25;
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    indd = 5+2*i-1; % for kmax = 3
    plot(dt:dt:N*dt, real(fuhat(indd,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_hat(indfmod(indd),1:N)), 'b', 'linewidth',2)
    title(['GB mode ( ', num2str(fkk(1,indd)),' , ', num2str(fkk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(fuhat(indd,:))+reshape(1*sqrt(sfu_post_cov(indcell,2*nqq+indd,:)),1,[]), real(fuhat(indd,end:-1:1))-reshape(1*sqrt(real(sfu_post_cov(indcell,2*nqq+indd,end:-1:1))),1,[])]','r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
    
    
    subplot(4,2,2*i)
    hold on
    indd = 24+2*i-1; % for kmax = 3
    plot(dt:dt:N*dt, real(fuhat(indd,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_hat(indfmod(indd),1:N)), 'b', 'linewidth',2)
    title(['GB mode ( ', num2str(fkk(1,indd)),' , ', num2str(fkk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(fuhat(indd,:))+reshape(1*sqrt(sfu_post_cov(indcell,2*nqq+indd,:)),1,[]), real(fuhat(indd,end:-1:1))-reshape(1*sqrt(real(sfu_post_cov(indcell,2*nqq+indd,end:-1:1))),1,[])]','r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')
end