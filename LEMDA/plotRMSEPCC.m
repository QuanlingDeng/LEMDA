%% plot rmse and pcc with respect to time in furier domain

file_name = sprintf('./uhat/ocn.mat', nqq);
load(file_name)
n = length(u_hat(:,1));

vec1 = real(u_hat(1,2000:end));

flerr = zeros(np, 1); flpcc = zeros(np,1);
for nqq = indv
    file_name = sprintf('./uhat/uhat%05d.mat', nqq);
    load(file_name)

    vec2 = real(u_post_mean(1,2000:end));

    rmse = sum( (vec1 - vec2).^2, 2);
    rrmse = sum( vec1.^2, 2);
    pcc = dot(vec1,vec2,2) ./ sqrt(dot(vec1,vec1,2) .* dot(vec2,vec2,2) );

    rmse = real(rmse); pcc = real(pcc); rrmse = real(rrmse);
    rrmse = sqrt(rmse./rrmse); rmse = sqrt(rmse);

    flerr(nqq) = mean(rrmse);
    flpcc(nqq) = mean(pcc);
end


figure
subplot(1,2,1)
hold on
plot(indv, flerr(indv), '-', 'linewidth',2)
title('Normalised RMSE','fontsize',14)
set(gca,'fontsize',24)
box on
xlabel('N_{obs}')


subplot(1,2,2)
hold on
plot(indv, flpcc(indv), '-', 'linewidth',2)
title('PCC','fontsize',14)
set(gca,'fontsize',24)
box on
xlabel('N_{obs}')

%%

figure
for i = 1:4
    subplot(4,2,2*i-1)
    hold on
    %indd = mod(24*(i-1)+1,40); % for kmax = 3
    indd = mod(60*(i-1)+1,100); % for kmax = 5
    plot(dt:dt:N*dt, u_hat(indd,1:N), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, u_post_mean(indd,:), 'r', 'linewidth',2)
    title(['(a) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(indd,:))+2*sqrt(real(u_post_cov(indd,:))), real(u_post_mean(indd,end:-1:1))-2*sqrt(real(u_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    set(gca,'fontsize',15)
    box on
    xlabel('t')


    subplot(4,2,2*i)
    hold on
    indd = Dim_Ug*2 + 15*(i-1)+4;
    plot(dt:dt:N*dt, u_hat(indd,1:N), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, u_post_mean(indd,:), 'r', 'linewidth',2)
    patch([dt:dt:N*dt,N*dt:-dt:dt], [real(u_post_mean(indd,:))+2*sqrt(real(u_post_cov(indd,:))), real(u_post_mean(indd,end:-1:1))-2*sqrt(real(u_post_cov(indd,end:-1:1)))],'r','facealpha',0.2,'linestyle','none')
    title(['(d) GB mode ( ', num2str(kk(1,indd)),' , ', num2str(kk(2,indd)), ' )'],'fontsize',14)

    set(gca,'fontsize',15)
    box on
    xlabel('t')
end





%% plot rmse and pcc with respect to time in physical domain

figure
for j = 0:3
    nqq = 2^j*10;
    file_name = sprintf('./err/err%05d.mat', nqq);
    load(file_name)

    subplot(2,2,1)
    hold on
    plot(dt:dt:N*dt, rrmse(:,1), '-', 'linewidth',2)
    title('Normalised RMSE','fontsize',14)
    set(gca,'fontsize',24)
    box on
    xlabel('t')

    subplot(2,2,3)
    hold on
    plot(dt:dt:N*dt, rrmse(:,2), '-', 'linewidth',2)
    title('Normalised RMSE','fontsize',14)
    set(gca,'fontsize',24)
    box on
    xlabel('t')


    subplot(2,2,2)
    hold on
    plot(dt:dt:N*dt, pcc(:,1), '-', 'linewidth',2)
    title('PCC','fontsize',14)
    set(gca,'fontsize',24)
    box on
    xlabel('t')

    subplot(2,2,4)
    hold on
    plot(dt:dt:N*dt, pcc(:,2), '-', 'linewidth',2)
    title('PCC','fontsize',14)
    set(gca,'fontsize',24)
    box on
    xlabel('t')

end

%% plot rmse and pcc with respect to number of observational particles in physical domain
indv = 10:10:200;
lerr = zeros(np, 2); lpcc = zeros(np,2);

for nqq = 10:10:200 %indv
    file_name = sprintf('./err/err%05d.mat', nqq);
    load(file_name)

    lerr(nqq, 1) = mean(rrmse(2001:end,1));
    lerr(nqq, 2) = mean(rrmse(2001:end,2));

    lpcc(nqq, 1) = mean(pcc(2001:end,1));
    lpcc(nqq, 2) = mean(pcc(2001:end,2));
end


figure
subplot(1,2,1)
hold on
plot(indv, lerr(indv,1), '-', 'linewidth',2)
plot(indv, lerr(indv,2), ':', 'linewidth',2)
title('Normalised RMSE','fontsize',14)
set(gca,'fontsize',24)
box on
xlabel('N_{obs}')


subplot(1,2,2)
hold on
plot(indv, lpcc(indv,1), '-', 'linewidth',2)
plot(indv, lpcc(indv,2), ':', 'linewidth',2)
title('PCC','fontsize',14)
set(gca,'fontsize',24)
box on
xlabel('N_{obs}')


return
np = 200;
for nqq = indv
    file_name = sprintf('./err2/err%05d.mat', nqq);
    load(file_name)
    
    rmse = real(rmse);
    rrmse = real(rrmse);
    pcc = real(pcc);

    rrmse = rrmse.*sqrt(rmse);
    save(['./err/err' num2str(nqq,'%05.f') '.mat'], "rmse", "rrmse", "pcc")
end
