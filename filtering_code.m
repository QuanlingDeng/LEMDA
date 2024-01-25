% Comparison of full and reduced order filter/smoother

% Lagrangian data assimilation with a given number of tracer L
rng(3); % fix the random number seed to reproduce results
sigma_xy = 0.01; % noise in the Lagrangian tracer equations
L = 20;80; % number of Lagrangian tracers
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
R0 = 0.01*eye(2*L+n,2*L+n); % initial value of posterior covariance (this cannot be zero if smoothing is later applied!)
u_post_mean = zeros(2*L+n,N); % posterior mean
u_post_mean(:,1) = mu0;
u_post_cov = zeros(2*L+n,N); % posterior covariance
u_post_cov(:,1) = diag(R0); % only save the diagonal elements

u_post_cov_all = zeros(2*L+n,2*L+n,N);
u_post_cov_all(:,:,1) = R0;
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
u_post_cov=u_post_cov/2;
R_end = R/2;
disp('Numerical value of the diagonal entry of R22 at final time')
disp(R(end,end))
R21 = R(2*L+1:end, 1:2*L); R12 = R21';
temp = diag(diag(R21*R12/sigma_xy^2));
% mu_s = zeros(2*L+n,N); 
% mu_s(:,end) = mu;
% R_s = zeros(2*L+n,N); % posterior covariance
% R_s(:,1) = diag(R); % only save the diagonal elements
% R_s_temp0 = R;
% for i = N-1:-1:1
%     Q = zeros(2*L,n);    
%     Q(1:L,:)     = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(L,1) * rk(1,:)));
%     Q(L+1:2*L,:) = beta* (exp(1i * x(:,i+1) * kk(1,:) + 1i * y(:,i+1) * kk(2,:)) .* (ones(L,1) * rk(2,:)));
%     a1_all = [-beta* eye(2*L,2*L), Q; zeros(n,2*L), a1];
%     mu_s(:,i) = mu_s(:,i+1) + (- a1_all*mu_s(:,i+1) + Sigma_u_all*Sigma_u_all' * (squeeze(u_post_cov_all(:,:,i+1)) \ (u_post_mean(:,i+1) - mu_s(:,i+1)))) * dt;
%     R_s_temp = R_s_temp0 +  ( - ( a1_all + (Sigma_u_all*Sigma_u_all') / (squeeze(u_post_cov_all(:,:,i+1))) ) * R_s_temp0 ...
%         - R_s_temp0 * ( a1_all' + (squeeze(u_post_cov_all(:,:,i+1)))\(Sigma_u_all*Sigma_u_all')  ) + Sigma_u_all*Sigma_u_all') * dt;
%     R_s_temp = (R_s_temp + R_s_temp')/2;
%     R_s(:,i) = diag(real(R_s_temp));
%     R_s_temp0 = R_s_temp;
% end
% pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Reduced filter %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reduced filter')
mu0 = [zeros(2*L,1);u_hat(:,1)]; % initial value of posterior mean
n = length(kk(1,:));
R0 = 0.001*eye(2*L+n,2*L+n); % initial value of posterior covariance (this cannot be zero if smoothing is later applied!)
u_post_mean_rd = zeros(2*L+n,N); % posterior mean
u_post_mean_rd(:,1) = mu0;
u_post_cov_rd = zeros(2*L+n,N); % posterior covariance
u_post_cov_rd(:,1) = diag(R0); % only save the diagonal elements
u_post_cov_all_rd = zeros(2*L+n,2*L+n,N);
u_post_cov_all_rd(:,:,1) = R0;
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
%     R22_temp = ((Sigma_u*Sigma_u')./(-a1 + sqrt(-a1 + L/sigma_xy^2*(Sigma_u*Sigma_u'))));
%     R22 = diag(diag(R22_temp));
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

 

 
figure
 
subplot(2,2,1)

hold on
plot(dt:dt:N*dt, real(u_hat(1,:)), 'b', 'linewidth',2)
plot(dt:dt:N*dt, real(u_post_mean(2*L+1,:)), 'r', 'linewidth',2)
plot(dt:dt:N*dt, real(u_post_mean_rd(2*L+1,:)), 'g', 'linewidth',2)
title(['(a) Filtering mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)        
post_upper = real(u_post_mean(2*L+1,:)) + 2 * sqrt(u_post_cov(2*L+1,:));
post_lower = real(u_post_mean(2*L+1,:)) - 2 * sqrt(u_post_cov(2*L+1,:));
post_upper_rd = real(u_post_mean_rd(2*L+1,:)) + 2 * sqrt(u_post_cov_rd(2*L+1,:));
post_lower_rd = real(u_post_mean_rd(2*L+1,:)) - 2 * sqrt(u_post_cov_rd(2*L+1,:));



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
% figure
% for i = 1:4
%     subplot(2,4,i)    
%     hold on
%     plot(dt:dt:N*dt, real(u(i,:)), 'b', 'linewidth',2)
%     plot(dt:dt:N*dt, real(u_post_mean(i,:)), 'r', 'linewidth',2)
%     plot(dt:dt:N*dt, real(u_post_mean_rd(i,:)), 'g', 'linewidth',2)
%     title(['Filtering mode ( ', num2str(kk(1,i)),' , ', num2str(kk(2,i)), ' )'],'fontsize',14)      
%     box on
%     set(gca,'fontsize',12)
%     if i == 1
%         legend('Truth','Full','Reduced')
%     end
%     xlim([1,9])
%     subplot(2,4,i+4)    
%     hold on
%     plot(dt:dt:N*dt, real(u(i,:)), 'b', 'linewidth',2)
%     plot(dt:dt:N*dt, real(mu_s(i,:)), 'r--', 'linewidth',2)
%     plot(dt:dt:N*dt, real(mu_s_rd(i,:)), 'g--', 'linewidth',2)
%     title(['Smoothing mode ( ', num2str(kk(1,i)),' , ', num2str(kk(2,i)), ' )'],'fontsize',14)      
%     box on
%     set(gca,'fontsize',12)
%     xlim([1,9])
% end

 
figure
for i = 1:3
    subplot(2,3,i)
    Dim_Grid = 50;
    [xx,yy] = meshgrid(linspace(-pi,pi,Dim_Grid), linspace(-pi,pi,Dim_Grid));
    x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
    u = exp(1i * x_vec * kk) * (u_hat(:,1+1000*(4*i-2)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+1000*(4*i-2)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    hold on
    contourf(xx, yy, sqrt(u.^2+ v.^2), 'linestyle','none')
    
    Dim_Grid = 25;
    [xx,yy] = meshgrid(linspace(-pi,pi,Dim_Grid), linspace(-pi,pi,Dim_Grid));
    x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
    u = exp(1i * x_vec * kk) * (u_hat(:,1+1000*(4*i-2)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+1000*(4*i-2)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    box on    
    title(['t = ', num2str(dt*1000*(4*i-2))])
    axis square
    plot(x(:,1+1000*(4*i-2)), y(:,1+1000*(4*i-2)),'ko','markersize',2,'linewidth',4);
    set(gca,'fontsize',12)
end
