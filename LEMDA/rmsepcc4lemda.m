%% RMSE root mean square error and PCC pattern cross correlation in physical domain
rmse = zeros(N,2); rrmse = zeros(N,2); pcc = zeros(N,2); 

% [xx,yy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx)); % 2nx or nx; no difference
% nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix
pnx = 2*ny; pny = 2*nx; fny = 2; fnx = 2; h = 2*pi/pny;
[xx,yy] = meshgrid(linspace(-pi+0.5*h,pi-0.5*h,pnx), linspace(-pi+0.5*h,pi-0.5*h,pny));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix

for j=2:N
    vx = exp(1i * nx_vec * kk) * (u_hat(:,j) .* transpose(rk(1,:)));
    vy = exp(1i * nx_vec * kk) * (u_hat(:,j) .* transpose(rk(2,:)));

    vxda = exp(1i * nx_vec * ckk) * (cu_post_mean(:,j) .* transpose(crk(1,:)));
    vyda = exp(1i * nx_vec * ckk) * (cu_post_mean(:,j) .* transpose(crk(2,:)));
    
%     fvx = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,j) .* transpose(frk(1,:)));
%     fvy = exp(1i * nx_vec * fkk) * (ffu_post_mean(2*nqq+1:end,j) .* transpose(frk(2,:)));
    
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
    vxda = real(vxda + reshape(fvx,[],1)); vyda = real(vyda + reshape(fvy,[],1)); 
    %vxda = real(vxda + fvx); vyda = real(vyda + fvy); 
    vx = real(vx); vy = real(vy); 
       
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
save(['./err/LEMDAnpi' num2str(npi,'%01.f') 'err.mat'], "err")

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