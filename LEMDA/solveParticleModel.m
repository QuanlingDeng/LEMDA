function maxo = solveParticleModel(domain, sig, sigv, np, dt, kk, rk, N, uhat, alpha)
% Solve the ice floe particle dynamics without collision force
% np: number of particles in domain
% nl: number of random particles we use for npagrangian DA data
% ny, ny, the mesh for Eulerian DA data (ny, nx) number of cells in domain
% sig: gives the noise strength
% alpha = 1; % ocean drag coefficients; 1/apha gives the decorrelation time as in a OU-process
npi = 1;
maxo = 0.0;

a = domain(2) - domain(1);
b = domain(4) - domain(3);

rng(2); % fix the random number seed to reproduce results
u = zeros(np,1);
v = zeros(np,1);
% x(:,1) = rand(np,1)*a - 0.5*a;
% y(:,1) = rand(np,1)*b - 0.5*a;
% % x(:,1) = 2*randn(np,1);
% % y(:,1) = 2*randn(np,1);
% x(:,1) = mod(real(x(:,1)) + 0.5*a, a) - 0.5*a; % periodic boundary conditions
% y(:,1) = mod(real(y(:,1)) + 0.5*b, b) - 0.5*a;

ns = ceil(sqrt(np));
h = 2*pi/ns;
[xx,yy] = meshgrid(linspace(-pi+0.5*h,pi-0.5*h,ns), linspace(-pi+0.5*h,pi-0.5*h,ns));
x(:,1) = reshape(xx,[],1);
y(:,1) = reshape(yy,[],1);
npp = randperm(ns^2);
x = x(npp(1:np));
y = y(npp(1:np));

savedn = 5000;
FloeX = zeros(np, savedn); FloeY = zeros(np, savedn); FloeU = zeros(np, savedn); FloeV = zeros(np, savedn); 

for i = 2:N+1 % generating the tracer locations
    x = x + u*dt + randn(np,1) * sig * sqrt(dt);
    y = y + v*dt + randn(np,1) * sig * sqrt(dt);

    uo = exp(1i * x * kk(1,:) + 1i * y * kk(2,:)) * (uhat(:,i-1) .* transpose(rk(1,:)));
    vo = exp(1i * x * kk(1,:) + 1i * y * kk(2,:)) * (uhat(:,i-1) .* transpose(rk(2,:)));

    u = u +  alpha * (uo - u) * dt + randn(np,1) * sigv * sqrt(dt);
    v = v +  alpha * (vo - v) * dt + randn(np,1) * sigv * sqrt(dt);
    
    x = mod(real(x) + 0.5*a, a) - 0.5*a; % periodic boundary conditions
    y = mod(real(y) + 0.5*b, b) - 0.5*a; % periodic boundary conditions
    
    modi = mod(i-2, savedn) + 1;
    FloeX(:,modi) = x; FloeY(:,modi) = y; FloeU(:,modi) = u; FloeV(:,modi) = v; 
    %La = [x y u v];
    if (mod(i, savedn)==1)
        %save(['/g/data/zv32/seaIceFloeData/cfdata/np' num2str(npi,'%02.f') 'time' num2str((i-1)/savedn,'%03.f') '.mat'],"FloeX", "FloeY","FloeU","FloeV");
        save(['./data/np' num2str(npi,'%02.f') 'time' num2str((i-1)/savedn,'%03.f') '.mat'],"FloeX", "FloeY","FloeU","FloeV");
    end
    
end


end