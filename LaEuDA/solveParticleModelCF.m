function maxo = solveParticleModelCF(domain, sig, sigv, np, dt, kk, rk, N, uhat, alpha)
% np: number of particles in domain
% nl: number of random particles we use for npagrangian DA data
% ny, ny, the mesh for Eulerian DA data (ny, nx) number of cells in domain
% sig: gives the noise strength
% alpha = 1; % ocean drag coefficients; 1/apha gives the decorrelation time as in a OU-process

maxo = 0.0;

a = domain(2) - domain(1);
b = domain(4) - domain(3);

rng(2); % fix the random number seed to reproduce results
r_min = 0.3*a / sqrt(np); r_max = r_min*2;  % assuming uniform size floes here
%thickness = ones(np,1); % thickness
r = r_min*ones(np,1); % radius
%r = r_min*(ones(np,1)+rand(np,1)); % radius

u = zeros(np,1);
v = zeros(np,1);
%x(:,1) = rand(np,1)*a - 0.5*a;
%y(:,1) = rand(np,1)*b - 0.5*a;
x(:,1) = 2*randn(np,1);
y(:,1) = 2*randn(np,1);
%x(:,1) = mod(real(x(:,1)) + 0.5*a, a) - 0.5*a; % periodic boundary conditions
%y(:,1) = mod(real(y(:,1)) + 0.5*b, b) - 0.5*a;

ns = 45;
h = 2*pi/ns;
[xx,yy] = meshgrid(linspace(-pi+0.5*h,pi-0.5*h,ns), linspace(-pi+0.5*h,pi-0.5*h,ns));
xt(:,1) = reshape(xx,[],1);
yt(:,1) = reshape(yy,[],1);
aa = randperm(2025); aa = sort(aa(1:np)); 
x(:,1) = xt(aa,1); y(:,1) = yt(aa,1); 

nx = ns; ny = nx; hx = a/nx; hy = b/ny; nc = 10;
cell4p = zeros(nx*ny, nc); % store the particle indices in each cell location; the first column gives the total number in this cell

for j=1:np
    indx = ceil((x(j) - domain(1))/hx);
    indy = ceil((y(j) - domain(3))/hy);
    ind = (indy-1)*nx + indx;
    
    cell4p(ind, 1) = cell4p(ind, 1) + 1;
    cell4p(ind, cell4p(ind, 1)+1) = j;
end

for i = 2:N+1 % generating the tracer locations
    
    CFx = zeros(np,1);
    CFy = zeros(np,1);
%     for j=1:np
%         floej = [r(j) x(j) y(j) u(j) v(j) 0.0];
%         for ki=j+1:np
%             floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
%             [CFtem, ~] = getCF(floej, floek);
%             CFx(j) = CFx(j) + CFtem(1);
%             CFy(j) = CFy(j) + CFtem(2); 
%             CFx(ki) = CFx(ki) - CFtem(1);
%             CFy(ki) = CFy(ki) - CFtem(2); 
%         end
%     end

    for j=1:np
        floej = [r(j) x(j) y(j) u(j) v(j) 0.0];
        
        indx = ceil((x(j) - domain(1))/hx);
        indy = ceil((y(j) - domain(3))/hy);
        ind = (indy-1)*nx + indx;
        
        for k=1:cell4p(ind, 1)
            ki = cell4p(ind, k+1);
            if ki~=j
                floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
                [CFtem, ~] = getCF(floej, floek);
                CFx(j) = CFx(j) + CFtem(1);
                CFy(j) = CFy(j) + CFtem(2);
            end
        end
        locCF;      
    end
    
    uo = exp(1i * x * kk(1,:) + 1i * y * kk(2,:)) * (uhat(:,i-1) .* transpose(rk(1,:)));
    vo = exp(1i * x * kk(1,:) + 1i * y * kk(2,:)) * (uhat(:,i-1) .* transpose(rk(2,:)));

    u = u +  alpha * (uo - u) * dt + CFx * dt  + randn(np,1) * sigv * sqrt(dt);
    v = v +  alpha * (vo - v) * dt + CFy * dt  + randn(np,1) * sigv * sqrt(dt);
    
    x = x + u*dt + randn(np,1) * sig * sqrt(dt);
    y = y + v*dt + randn(np,1) * sig * sqrt(dt);
    x = mod(real(x) + 0.5*a, a) - 0.5*a; % periodic boundary conditions
    y = mod(real(y) + 0.5*b, b) - 0.5*a; % periodic boundary conditions
    
    cell4p = 0.0*cell4p;
    for j=1:np
        indx = ceil((x(j) - domain(1))/hx);
        indy = ceil((y(j) - domain(3))/hy);
        ind = (indy-1)*nx + indx;
        
        cell4p(ind, 1) = cell4p(ind, 1) + 1;
        cell4p(ind, cell4p(ind, 1)+1) = j;
    end
    
    La = [x y u v];
    save(['./data/time' num2str(i-1,'%05.f') '.mat'],"La");

%     mm = max([abs(u); abs(v)]);
%     if maxo< mm
%         maxo = mm;
%     end
end

end