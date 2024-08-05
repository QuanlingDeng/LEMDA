% Code: Generating the true flow field for DA

% 2D incompressible flow with random amplitude based on linear stochastic
% models, mimicking the QG flow.
K_max = 4; % the range of Fourier modes is [-K_max, K_max]^2
k = zeros(2, (2 * K_max + 1) * (2 * K_max + 1)); % Total number of Fourier wavenumbers

% arranging Fourier wavenumbers
% arranging in such a way that the complex conjugates modes are next to
% each other, namely (-k1,-k2) will be next to (k1,k2). This will
% facilitate data assimilation and allows it to be a nearly block diagonal
% matrix where each block is 2 by 2.
m = 1;
for i = - K_max : K_max
    if i < 0
        for j = - K_max: i
            k(1, m) = i;
            k(2, m) = j;
            m = m + 2;
        end
    else
        for j = - K_max : i - 1
            k(1, m) = i;
            k(2, m) = j;
            m = m + 2;
        end
    end
end
k(:, 2: 2: end - 1) = - k(:, 1 : 2 : end - 2);

% showing the grids of Fourier wavenumbers
figure
hold on
plot(k(1,1:2:end-2), k(2,1:2:end-2), 'ro', 'linewidth',4);
plot(k(1,2:2:end-1), k(2,2:2:end-1), 'go', 'linewidth',4);
plot(k(1,end), k(2,end), 'ko', 'linewidth',4);
box on
grid on
set(gca, 'fontsize', 12)
title('Fourier wavenumbers', 'fontsize', 14)
xlabel('k_1')
ylabel('k_2')

% wavenumbers
kk = k(:,1:end-1); % remove the (0,0) mode by assuming no background flow.

% eigenvectors 
rk1 = [1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (-1i * k(2,:));
    1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (1i * k(1,:))];
rk1 = rk1(:,1:end-1);
rk = rk1;


% stochastic systems for each Fourier mode
% written in the vector form
% dU = (a1 U + a0 ) dt + b1 dW

rng(100) % fix the random number seed
T = 10.0001; % total time
dt = 0.0001; % time step
N = round(T/dt); % total number of time steps within the given time interval
Dim_U = length(kk(1,:)); % dimension of the system
Dim_UB = Dim_U;
u_hat = zeros(Dim_U,N); % define all the Fourier modes
d_B = 1; % damping  
sigma_B = 0.05; % noise 


% b1: noise coefficient; a1: damping and phase 
b1 = zeros(Dim_U, Dim_U);
a1 = - eye(Dim_UB) * d_B;  
L_u = a1;
for i = 1: 2: Dim_UB
    b1(i,i) = 1 / sqrt(2) * sigma_B;
    b1(i+1, i+1) = -1i / sqrt(2) * sigma_B;
    b1(i, i+1) = 1i / sqrt(2) * sigma_B;
    b1(i+1, i) = 1 / sqrt(2) * sigma_B;
end
Sigma_u = b1;

% numerical integration
rd = randn(Dim_U,N);

a0 = zeros(Dim_U, 1);
for i = 2:N
    t = i*dt; 
    u_hat(:,i) = u_hat(:,i-1) + (L_u * u_hat(:,i-1) + a0) * dt + Sigma_u * sqrt(dt) * rd(:, i);
end
 

% spatiotemporal reconstruction
Dim_Grid = 25;
[xx,yy] = meshgrid(linspace(-pi,pi,Dim_Grid), linspace(-pi,pi,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
nn = round(T/dt/100)/2;
uu_all = zeros(1,nn);

% plotting the spatial flow field
figure
for i = 1:9
    subplot(3,3,i)
    u = exp(1i * x_vec * kk) * (u_hat(:,1000*i) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1000*i) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    hold on
    contourf(xx,yy, sqrt(u.^2+v.^2),20);
    colormap jet
    shading interp
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([-pi, pi ])
    ylim([-pi, pi ])
    box on    
    title(['t = ', num2str(dt*1000*i)])
    colorbar
    axis square
end

% plotting the time series of different Fourier modes
figure
for i = 1:4
    subplot(2,2,i)
    if i == 1
        plot(dt:dt:N*dt, u_hat(1,:), 'b', 'linewidth',2)
        title(['(a) mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 2
        plot(dt:dt:N*dt, u_hat(2,:), 'b', 'linewidth',2)
        title(['(b) mode ( ', num2str(kk(1,2)),' , ', num2str(kk(2,2)), ' )'],'fontsize',14)
    elseif i == 3
        plot(dt:dt:N*dt, u_hat(3,:), 'b', 'linewidth',2)
        title(['(c) mode ( ', num2str(kk(1,3)),' , ', num2str(kk(2,3)), ' )'],'fontsize',14)
    elseif i == 4
        plot(dt:dt:N*dt, u_hat(4,:), 'b', 'linewidth',2)
        title(['(d) mode ( ', num2str(kk(1,4)),' , ', num2str(kk(2,4)), ' )'],'fontsize',14)
    end
    set(gca,'fontsize',12)
    box on
    xlabel('t')
end