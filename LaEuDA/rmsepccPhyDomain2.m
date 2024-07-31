%% RMSE root mean square error and PCC pattern cross correlation in physical domain
rmse = zeros(N,2); rrmse = zeros(N,2); pcc = zeros(N,2); 

[xx,yy] = meshgrid(linspace(-pi+0.5*hx,pi-0.5*hx,nx), linspace(-pi+0.5*hx,pi-0.5*hx,nx)); % 2nx or nx; no difference
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)]; % becoming a two column matrix

for j=1:N
%     vx = exp(1i * nx_vec * kk) * (u_hat(:,j) .* transpose(rk(1,:)));
%     vy = exp(1i * nx_vec * kk) * (u_hat(:,j) .* transpose(rk(2,:)));
% 
%     vxda = exp(1i * nx_vec * kk) * (u_post_mean(:,j) .* transpose(rk(1,:)));
%     vyda = exp(1i * nx_vec * kk) * (u_post_mean(:,j) .* transpose(rk(2,:)));
%     vxda = exp(1i * nx_vec * kk) * (u_post_mean(1:Dim_U,j) .* transpose(rk(1,:))); 
%     vyda = exp(1i * nx_vec * kk) * (u_post_mean(1:Dim_U,j) .* transpose(rk(2,:)));
%     vxda = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,j) .* transpose(rk(1,:))); 
%     vyda = exp(1i * nx_vec * kk) * (u_post_mean(2*nqq+1:end,j) .* transpose(rk(2,:)));

%     vx = real(vx); vy = real(vy); 
%     vxda = real(vxda); vyda = real(vyda); 
    
    vx = real(u_hat(:,j)); vy = imag(u_hat(:,j)); 
    %vxda = real(u_post_mean(:,j)); vyda = imag(u_post_mean(:,j)); 
    vxda = real(u_post_mean(2*nqq+1:end,j)); vyda = imag(u_post_mean(2*nqq+1:end,j)); 

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

%%
figure
% subplot(1,2,1)
% hold on
% plot(dt:dt:N*dt, rmse(:,1), '-b', 'linewidth',2)
% plot(dt:dt:N*dt, rmse(:,2), '--r', 'linewidth',2)
% %plot(dt:dt:N*dt, rrmse(:,1), '-r', 'linewidth',2)
% %plot(dt:dt:N*dt, rrmse(:,2), '--r', 'linewidth',2)
% title('RMSE','fontsize',14)
% set(gca,'fontsize',16)
% legend('vx','vy')
% box on
% xlabel('t')

subplot(1,2,1)
hold on
%plot(dt:dt:N*dt, rmse(:,1), '-b', 'linewidth',2)
%plot(dt:dt:N*dt, rmse(:,2), '--b', 'linewidth',2)
% plot(dt:dt:N*dt, rrmse(:,1), '-b', 'linewidth',2)
% plot(dt:dt:N*dt, rrmse(:,2), '--r', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*rrmse(:,1).^2+0.5*rrmse(:,2).^2), '-b', 'linewidth',2)
title('Normalised RMSE','fontsize',14)
set(gca,'fontsize',16)
%legend('vx','vy')
box on
xlabel('t')

subplot(1,2,2)
hold on
% plot(dt:dt:N*dt, pcc(:,1), '-b', 'linewidth',2)
% plot(dt:dt:N*dt, pcc(:,2), '--r', 'linewidth',2)
plot(dt:dt:N*dt, sqrt(0.5*pcc(:,1).^2+0.5*pcc(:,2).^2), '-b', 'linewidth',2)
title('PCC','fontsize',14)
set(gca,'fontsize',16)
%legend('vx','vy')
box on
xlabel('t')