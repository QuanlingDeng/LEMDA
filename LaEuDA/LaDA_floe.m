% Code: Lagrangian data assimilation with a given number of tracer L
% Note that the code also generates the smoother solution although it was
% not included in the paper

Generating_GB_flow % Generating ocean current

rng(300); % fix the random number seed to reproduce results
sigma_xy = 0.001; % noise in the Lagrangian tracer equations
L = 10; % number of Lagrangian tracers
x = zeros(L,N); % tracer displacement
y = zeros(L,N);
u = zeros(L,N); % tracer velocity
v = zeros(L,N);
x(:,1) = rand(L,1)*2*pi;
y(:,1) = rand(L,1)*2*pi;
beta = 1; % the drag coefficient
for i = 2:N % generating the tracer locations
    x(:,i) = x(:,i-1) + u(:,i-1) * dt + randn(L,1) * sigma_xy * sqrt(dt); 
    u(:,i) = u(:,i-1) + beta * (real(exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) * (u_hat(:,i-1) .* transpose(rk(1,:)))) - u(:,i-1)) * dt;
    y(:,i) = y(:,i-1) + v(:,i-1) * dt + randn(L,1) * sigma_xy * sqrt(dt); 
    v(:,i) = v(:,i-1) + beta * (real(exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) * (u_hat(:,i-1) .* transpose(rk(2,:)))) - v(:,i-1)) * dt;
    x(:,i) = mod(real(x(:,i)) + pi, 2*pi) - pi; % periodic boundary conditions
    y(:,i) = mod(real(y(:,i)) + pi, 2*pi) - pi; % periodic boundary conditions
end

l = length(k(1,:)); % number of Fourier wavenumbers
InvBoB = eye(2*L)/sigma_xy/sigma_xy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Full filter %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Full filter')
mu0 = [zeros(2*L,1);u_hat(:,1)]; % initial value of posterior mean
n = length(kk(1,:));
R0 = 0.001*eye(2*L+n,2*L+n); % initial value of posterior covariance (this cannot be zero if smoothing is later applied!)
u_post_mean = zeros(2*L+n,N); % posterior mean
u_post_mean(:,1) = mu0;
u_post_cov = zeros(2*L+n,N); % posterior covariance
u_post_cov(:,1) = diag(R0); % only save the diagonal elements

u_post_cov_all = zeros(2*L+n,2*L+n,N);
u_post_cov_all(:,:,1) = R0;
% solving the filtering solution
for i = 2:N
    x0 = x(:,i-1); 
    y0 = y(:,i-1); 
    x1 = x(:,i); 
    y1 = y(:,i); 
    x_diff = x1-x0;
    y_diff = y1-y0;
    % need to take into account the periodic boundary conditions
    x_diff(x_diff > pi)  = x_diff(x_diff>pi)  - 2 * pi;
    x_diff(x_diff < -pi) = x_diff(x_diff<-pi) + 2 * pi;
    y_diff(y_diff > pi)  = y_diff(y_diff>pi)  - 2 * pi;
    y_diff(y_diff < -pi) = y_diff(y_diff<-pi) + 2 * pi;
    % matrix for filtering
    A1 = [eye(2*L), zeros(2*L, n)];
    Q = zeros(2*L,n);    
    Q(1:L,:)     = beta* (exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(L,1) * rk(1,:)));
    Q(L+1:2*L,:) = beta* (exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(L,1) * rk(2,:)));
    a1_all = [-beta* eye(2*L,2*L), Q; zeros(n,2*L), a1];
    a0_all = [zeros(2*L,1);a0];
    Sigma_u_all = zeros(2*L+n,2*L+n);
    Sigma_u_all(2*L+1:end, 2*L+1:end) = Sigma_u;
    % update the posterior mean and posterior covariance
    mu = mu0 + (a0_all + a1_all * mu0) * dt + (R0 * A1') * InvBoB * ([x_diff;y_diff] - A1 * mu0 * dt);
    R =  R0 + (a1_all * R0 + R0* a1_all' + Sigma_u_all*Sigma_u_all' - (R0*A1') * InvBoB * (R0*A1')')*dt;
    R = (R+R')/2;
    u_post_mean(:,i) = mu;
    u_post_cov(:,i) = diag(real(R));
    mu0 = mu;
    R0 = R;
    u_post_cov_all(:,:,i) = R;
end
R_end = R;
% solving the smoothing solution 
R21 = R(2*L+1:end, 1:2*L); R12 = R21';
temp = diag(diag(R21*R12/sigma_xy^2));
mu_s = zeros(2*L+n,N); 
mu_s(:,end) = mu;
R_s = zeros(2*L+n,N); % posterior covariance
R_s(:,1) = diag(R); % only save the diagonal elements
R_s_temp0 = R;
for i = N-1:-1:1
    Q = zeros(2*L,n);    
    Q(1:L,:)     = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(L,1) * rk(1,:)));
    Q(L+1:2*L,:) = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(L,1) * rk(2,:)));
    a1_all = [-beta* eye(2*L,2*L), Q; zeros(n,2*L), a1];
    mu_s(:,i) = mu_s(:,i+1) + (- a1_all*mu_s(:,i+1) + Sigma_u_all*Sigma_u_all' * (squeeze(u_post_cov_all(:,:,i+1)) \ (u_post_mean(:,i+1) - mu_s(:,i+1)))) * dt;
    R_s_temp = R_s_temp0 +  ( - ( a1_all + (Sigma_u_all*Sigma_u_all') / (squeeze(u_post_cov_all(:,:,i+1))) ) * R_s_temp0 ...
        - R_s_temp0 * ( a1_all' + (squeeze(u_post_cov_all(:,:,i+1)))\(Sigma_u_all*Sigma_u_all')  ) + Sigma_u_all*Sigma_u_all') * dt;
    R_s_temp = (R_s_temp + R_s_temp')/2;
    R_s(:,i) = diag(real(R_s_temp));
    R_s_temp0 = R_s_temp;
end

%% computing the skill scores
Dim_Grid = 25;
[xx,yy] = meshgrid(linspace(-pi,pi,Dim_Grid), linspace(-pi,pi,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
nn = round(T/dt/100)/2;
uu_all = zeros(1,nn);

Relative_Entropy = zeros(1,round(N/10));
Dim = length(u_post_mean(2*L+1:end,end));
R_eq = sigma_B.^2/2/d_B * eye(length(u_post_mean(2*L+1:end,end)));
mu_eq = zeros(size(u_post_mean(2*L+1:end,end)));
RMSE_temp = zeros(1,round(N/10));
Corr_temp = zeros(1,round(N/10));
for i = 1:N/10
    u_phy_truth = real(exp(1i * x_vec * kk) * (u_hat(:,i*10) .* transpose(rk(1,:))));
    u_phy_post  = real(exp(1i * x_vec * kk) * (u_post_mean(2*L+1:end,i*10) .* transpose(rk(1,:))));
    RMSE_temp(i) = sqrt(sum((u_phy_truth - u_phy_post).^2)/Dim_Grid^2)/std(u_phy_truth);
    temp = corrcoef(u_phy_truth,u_phy_post);
    Corr_temp(i) = temp(1,2);
    Relative_Entropy(i) = real(1/2 * ( (u_post_mean(2*L+1:end,i*10)-mu_eq)' / R_eq * (u_post_mean(2*L+1:end,i*10)-mu_eq) )) + ...
                         real(1/2 * ( sum((u_post_cov(2*L+1:end,i*10))./diag(R_eq)) - Dim - log(det( diag(u_post_cov(2*L+1:end,i*10))*inv((R_eq)))) ));
end
RMSE_full = mean(RMSE_temp);
Corr_full = mean(Corr_temp);
Info_full = mean(Relative_Entropy);
disp(['RMSE, Corr and Uncertainty Reduction in full filter:'])
disp([RMSE_full, Corr_full, Info_full])
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Reduced filter %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reduced filter')
search_r_equil
mu0 = [zeros(2*L,1);u_hat(:,1)]; % initial value of posterior mean
n = length(kk(1,:));
R0 = 0.01*eye(2*L+n,2*L+n); % initial value of posterior covariance (this cannot be zero if smoothing is later applied!)
u_post_mean_rd = zeros(2*L+n,N); % posterior mean
u_post_mean_rd(:,1) = mu0;
u_post_cov_rd = zeros(2*L+n,N); % posterior covariance
u_post_cov_rd(:,1) = diag(R0); % only save the diagonal elements
u_post_cov_all_rd = zeros(2*L+n,2*L+n,N);
u_post_cov_all_rd(:,:,1) = R0;
% solving the filtering solution
for i = 2:N
    x0 = x(:,i-1); 
    y0 = y(:,i-1); 
    x1 = x(:,i); 
    y1 = y(:,i); 
    x_diff = x1-x0;
    y_diff = y1-y0;
    % need to take into account the periodic boundary conditions
    x_diff(x_diff > pi)  = x_diff(x_diff>pi)  - 2 * pi;
    x_diff(x_diff < -pi) = x_diff(x_diff<-pi) + 2 * pi;
    y_diff(y_diff > pi)  = y_diff(y_diff>pi)  - 2 * pi;
    y_diff(y_diff < -pi) = y_diff(y_diff<-pi) + 2 * pi;
    % matrix for filtering
    A1 = [eye(2*L), zeros(2*L, n)];
    Q = zeros(2*L,n);    
    Q(1:L,:)     = beta* (exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(L,1) * rk(1,:)));
    Q(L+1:2*L,:) = beta* (exp(1i * x(:,i-1) * kk(1,:) + 1i * y(:,i-1) * kk(2,:)) .* (ones(L,1) * rk(2,:)));
    a1_all = [-beta* eye(2*L,2*L), Q; zeros(n,2*L), a1];
    a0_all = [zeros(2*L,1);a0];
    Sigma_u_all = zeros(2*L+n,2*L+n);
    Sigma_u_all(2*L+1:end, 2*L+1:end) = Sigma_u;
    % update the posterior mean and posterior covariance
    R11_0 = R0(1:2*L, 1:2*L);
    R21_0 = R0(2*L+1:end, 1:2*L); R12_0 = R21_0';
    R22_0 = R0(2*L+1:end, 2*L+1:end);
    R22 = r2_rd*eye(n);%(- temp + Sigma_u * Sigma_u')/2 * inv(-a1);
    R11 = R11_0 + (-2*beta * R11_0 + Q * R21_0 + R12_0 * Q' - 1/sigma_xy^2 * R11_0 * R11_0) * dt;
    R12 = R12_0 + (-beta*R12_0 + Q* R22_0- R12_0* (-a1)- 1/sigma_xy^2 * R11_0 * R12_0) * dt;
    R21 = R12';
    R = [R11,R12;R21,R22];
%     R =  R0 + (a1_all * R0 + R0* a1_all' + Sigma_u_all*Sigma_u_all' - (R0*A1') * InvBoB * (R0*A1')')*dt;
%     R(1:2*L,1:2*L) = R(1:2*L,1:2*L) .* eye(2*L);
%     R(2*L+1:end, 2*L+1:end) = R(2*L+1:end, 2*L+1:end).*eye(n);
    mu = mu0 + (a0_all + a1_all * mu0) * dt + (R0 * A1') * InvBoB * ([x_diff;y_diff] - A1 * mu0 * dt);
    
    u_post_mean_rd(:,i) = mu;
    u_post_cov_rd(:,i) = diag(real(R));
    mu0 = mu;
    R0 = R;
    u_post_cov_all_rd(:,:,i) = R;
end
R_end_rd = R;
% solving the smoothing solution
mu_s_rd = zeros(2*L+n,N); 
mu_s_rd(:,end) = mu;
R_s_rd = zeros(2*L+n,N); % posterior covariance
R_s_rd(:,1) = diag(R); % only save the diagonal elements
R_s_temp0 = R;
for i = N-1:-1:1
    Q = zeros(2*L,n);    
    Q(1:L,:)     = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(L,1) * rk(1,:)));
    Q(L+1:2*L,:) = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(L,1) * rk(2,:)));
    a1_all = [-beta* eye(2*L,2*L), Q; zeros(n,2*L), a1];
    mu_s_rd(:,i) = mu_s_rd(:,i+1) + (- a1_all*mu_s_rd(:,i+1) + Sigma_u_all*Sigma_u_all' * (squeeze(u_post_cov_all_rd(:,:,i+1)) \ (u_post_mean_rd(:,i+1) - mu_s_rd(:,i+1)))) * dt;
    R_s_temp = R_s_temp0 +  ( - ( a1_all + (Sigma_u_all*Sigma_u_all') / (squeeze(u_post_cov_all_rd(:,:,i+1))) ) * R_s_temp0 ...
        - R_s_temp0 * ( a1_all' + (squeeze(u_post_cov_all(:,:,i+1)))\(Sigma_u_all*Sigma_u_all') ) + Sigma_u_all*Sigma_u_all') * dt;
    R_s_temp = (R_s_temp + R_s_temp')/2;
    R_s_rd(:,i) = diag(real(R_s_temp));
    R_s_temp0 = R_s_temp;
end
% computing the skill scores
Relative_Entropy = zeros(1,round(N/10));
Dim = length(u_post_mean_rd(2*L+1:end,end));
R_eq = sigma_B.^2/2/d_B * eye(length(u_post_mean_rd(2*L+1:end,end)));
mu_eq = zeros(size(u_post_mean_rd(2*L+1:end,end)));
RMSE_temp = zeros(1,round(N/10));
Corr_temp = zeros(1,round(N/10));
for i = 1:N/10
    u_phy_truth = real(exp(1i * x_vec * kk) * (u_hat(:,i*10) .* transpose(rk(1,:))));
    u_phy_post_rd  = real(exp(1i * x_vec * kk) * (u_post_mean_rd(2*L+1:end,i*10) .* transpose(rk(1,:))));
    RMSE_temp(i) = sqrt(sum((u_phy_truth - u_phy_post_rd).^2)/Dim_Grid^2)/std(u_phy_truth);
    temp = corrcoef(u_phy_truth,u_phy_post_rd);
    Corr_temp(i) = temp(1,2);
    Relative_Entropy(i) = real(1/2 * ( (u_post_mean_rd(2*L+1:end,i*10)-mu_eq)' / R_eq * (u_post_mean_rd(2*L+1:end,i*10)-mu_eq) )) + ...
                         real(1/2 * ( sum((u_post_cov_rd(2*L+1:end,i*10))./diag(R_eq)) - Dim - log(det( diag(u_post_cov_rd(2*L+1:end,i*10))*inv((R_eq)))) ));
end
RMSE_rd = mean(RMSE_temp);
Corr_rd = mean(Corr_temp);
Info_rd = mean(Relative_Entropy);
disp(['RMSE, Corr and Uncertainty Reduction in reduced filter:'])
disp([RMSE_rd, Corr_rd, Info_rd])

% showing the time series of different modes 
figure
for i = [1,3]
    subplot(2,2,i)
    if i == 1
        hold on
        plot(dt:dt:N*dt, real(u_hat(1,:)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean(2*L+1,:)), 'r', 'linewidth',2)
        plot(dt:dt:N*dt, real(u_post_mean_rd(2*L+1,:)), 'g', 'linewidth',2)
        title(['(a) Filtering mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)        
        post_upper = real(u_post_mean(2*L+1,:)) + 2 * sqrt(u_post_cov(2*L+1,:));
        post_lower = real(u_post_mean(2*L+1,:)) - 2 * sqrt(u_post_cov(2*L+1,:));
        post_upper_rd = real(u_post_mean_rd(2*L+1,:)) + 2 * sqrt(u_post_cov_rd(2*L+1,:));
        post_lower_rd = real(u_post_mean_rd(2*L+1,:)) - 2 * sqrt(u_post_cov_rd(2*L+1,:));

    elseif i == 3
        hold on
        plot(dt:dt:N*dt, real(u_hat(1,:)), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, real(mu_s(2*L+1,:)), 'r', 'linewidth',2)
        plot(dt:dt:N*dt, real(mu_s_rd(2*L+1,:)), 'g', 'linewidth',2)
        title(['(b) Smoothing mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)        
        post_upper = real(mu_s(2*L+1,:)) + 2 * sqrt(R_s(2*L+1,:));
        post_lower = real(mu_s(2*L+1,:)) - 2 * sqrt(R_s(2*L+1,:));
        post_upper_rd = real(mu_s_rd(2*L+1,:)) + 2 * sqrt(R_s_rd(2*L+1,:));
        post_lower_rd = real(mu_s_rd(2*L+1,:)) - 2 * sqrt(R_s_rd(2*L+1,:));

    end
    
    tt = dt:dt:N*dt;
    patch([tt,tt(end:-1:1)],[post_lower,post_upper(end:-1:1)],'r','facealpha',0.2,'linestyle','none');
    patch([tt,tt(end:-1:1)],[post_lower_rd,post_upper_rd(end:-1:1)],'g','facealpha',0.2,'linestyle','none');
    set(gca,'fontsize',12)
    if i == 1
        legend('Truth','Full','Reduced','2std Full','2std Reduced')
        yylim = get(gca,'ylim');
    end
    box on
    xlabel('t')
    xlim([1,9])
    ylim([yylim(1),yylim(2)])
end
% showing the covariance in the full filter 
subplot(2,2,2)
pcolor(real(R_end))
colorbar
title('(c) Full cov at t = 10 (filter)')
box on
set(gca,'fontsize',12)
axis equal
c = get(gca,'clim');
hold on 
plot([2*L,2*L],[2*L,2*L+n],'r')
plot([2*L,2*L+n],[2*L,2*L],'r')
% showing the covariance in the reduced filter 
subplot(2,2,4)
pcolor(real(R_end_rd))
colorbar
title('(d) Reduced-order cov at t = 10 (filter)')
box on
set(gca,'fontsize',12)
axis equal
caxis([c(1),c(2)])
hold on 
plot([2*L,2*L],[2*L,2*L+n],'r')
plot([2*L,2*L+n],[2*L,2*L],'r')

% plotting the tracer velocity field (not shown in the paper)
figure
for i = 1:4
    subplot(2,4,i)    
    hold on
    plot(dt:dt:N*dt, real(u(i,:)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean(i,:)), 'r', 'linewidth',2)
    plot(dt:dt:N*dt, real(u_post_mean_rd(i,:)), 'g', 'linewidth',2)
    title(['Filtering mode ( ', num2str(kk(1,i)),' , ', num2str(kk(2,i)), ' )'],'fontsize',14)      
    box on
    set(gca,'fontsize',12)
    if i == 1
        legend('Truth','Full','Reduced')
    end
    xlim([1,9])
    subplot(2,4,i+4)    
    hold on
    plot(dt:dt:N*dt, real(u(i,:)), 'b', 'linewidth',2)
    plot(dt:dt:N*dt, real(mu_s(i,:)), 'r--', 'linewidth',2)
    plot(dt:dt:N*dt, real(mu_s_rd(i,:)), 'g--', 'linewidth',2)
    title(['Smoothing mode ( ', num2str(kk(1,i)),' , ', num2str(kk(2,i)), ' )'],'fontsize',14)      
    box on
    set(gca,'fontsize',12)
    xlim([1,9])
end

 


% Comparing the flow fields
figure 
subplot(2,3,1)
u = exp(1i * x_vec * kk) * (u_hat(:,8000) .* transpose(rk(1,:)));
v = exp(1i * x_vec * kk) * (u_hat(:,8000) .* transpose(rk(2,:)));
u = reshape(u, Dim_Grid, Dim_Grid);
v = reshape(v, Dim_Grid, Dim_Grid);
hold on
contourf(xx,yy, sqrt(u.^2+v.^2),20);
colormap jet
shading interp
quiver(xx, yy, u, v, 'linewidth',1)
xlim([0, 2*pi ])
ylim([0, 2*pi ])
xlim([-pi, pi ])
ylim([-pi, pi ])
box on    
title(['Truth; t = ', num2str(8000*dt)])
set(gca,'fontsize',12) 

subplot(2,3,2)
u = real(exp(1i * x_vec * kk) * (u_post_mean(2*L+1:end,8000) .* transpose(rk(1,:))));
v = real(exp(1i * x_vec * kk) * (u_post_mean(2*L+1:end,8000) .* transpose(rk(2,:))));
u = reshape(u, Dim_Grid, Dim_Grid);
v = reshape(v, Dim_Grid, Dim_Grid);
hold on
contourf(xx,yy, sqrt(u.^2+v.^2),20);
colormap jet
shading interp
quiver(xx, yy, u, v, 'linewidth',1)
xlim([0, 2*pi ])
ylim([0, 2*pi ])
xlim([-pi, pi ])
ylim([-pi, pi ])
box on    
title(['Full filter; t = ', num2str(8000*dt)])
set(gca,'fontsize',12) 

subplot(2,3,3)
u = real(exp(1i * x_vec * kk) * (u_post_mean_rd(2*L+1:end,8000) .* transpose(rk(1,:))));
v = real(exp(1i * x_vec * kk) * (u_post_mean_rd(2*L+1:end,8000) .* transpose(rk(2,:))));
u = reshape(u, Dim_Grid, Dim_Grid);
v = reshape(v, Dim_Grid, Dim_Grid);
hold on
contourf(xx,yy, sqrt(u.^2+v.^2),20);
colormap jet
shading interp
quiver(xx, yy, u, v, 'linewidth',1)
xlim([0, 2*pi ])
ylim([0, 2*pi ])
xlim([-pi, pi ])
ylim([-pi, pi ])
box on    
title(['Reduced filter; t = ', num2str(8000*dt)])
set(gca,'fontsize',12) 

subplot(2,3,5)
u = real(exp(1i * x_vec * kk) * (mu_s(2*L+1:end,8000) .* transpose(rk(1,:))));
v = real(exp(1i * x_vec * kk) * (mu_s(2*L+1:end,8000) .* transpose(rk(2,:))));
u = reshape(u, Dim_Grid, Dim_Grid);
v = reshape(v, Dim_Grid, Dim_Grid);
hold on
contourf(xx,yy, sqrt(u.^2+v.^2),20);
colormap jet
shading interp
quiver(xx, yy, u, v, 'linewidth',1)
xlim([0, 2*pi ])
ylim([0, 2*pi ])
xlim([-pi, pi ])
ylim([-pi, pi ])
box on    
title(['Full smoother; t = ', num2str(8000*dt)])
set(gca,'fontsize',12) 

subplot(2,3,6)
u = real(exp(1i * x_vec * kk) * (mu_s_rd(2*L+1:end,8000) .* transpose(rk(1,:))));
v = real(exp(1i * x_vec * kk) * (mu_s_rd(2*L+1:end,8000) .* transpose(rk(2,:))));
u = reshape(u, Dim_Grid, Dim_Grid);
v = reshape(v, Dim_Grid, Dim_Grid);
hold on
contourf(xx,yy, sqrt(u.^2+v.^2),20);
colormap jet
shading interp
quiver(xx, yy, u, v, 'linewidth',1)
xlim([0, 2*pi ])
ylim([0, 2*pi ])
xlim([-pi, pi ])
ylim([-pi, pi ])
box on    
title(['Reduced smoother; t = ', num2str(8000*dt)])
set(gca,'fontsize',12) 
