% generate noise coefficient for the EuDA sigma_n = 0.15; % noise in the number density discretized ODEs
function sig = getNoise4EuDA(domain, np, fnx, nx, dt, T)

n2 = floor(fnx/nx); % size ratio
Nt = floor(T/dt);
beta = 1.0;
ns = Nt;

file_name = sprintf('./uhat/ocn.mat');
load(file_name)

% coarse scale 
ny = nx;
hx = (domain(2) - domain(1))/nx;
hy = (domain(4) - domain(3))/ny;
% npsol0 = zeros(ny, nx);
% npusol0 = zeros(ny, nx);
% npvsol0 = zeros(ny, nx);
npsol = zeros(ny, nx);
npusol = zeros(ny, nx);
npvsol = zeros(ny, nx);

% fine scale 
fny = fnx;
fhx = (domain(2) - domain(1))/fnx;
fhy = (domain(4) - domain(3))/fny;
% fnpsol0 = zeros(fny, fnx);
% fnpusol0 = zeros(fny, fnx);
% fnpvsol0 = zeros(fny, fnx);
fnpsol = zeros(fny, fnx);
fnpusol = zeros(fny, fnx);
fnpvsol = zeros(fny, fnx);

% file_name = sprintf('./data/time%05d.mat', Nt-1);
% load(file_name)
% 
% indx = ceil( ( La(:,1) - domain(1) )/hx );
% indy = ceil( ( La(:,2)  - domain(3) )/hy );
% findx = ceil( ( La(:,1) - domain(1) )/fhx );
% findy = ceil( ( La(:,2)  - domain(3) )/fhy );
% 
% for k=1:np
%     npsol0(indy(k), indx(k)) = npsol0(indy(k), indx(k)) + 1;
%     npusol0(indy(k), indx(k)) = npusol0(indy(k), indx(k)) + La(k,3);
%     npvsol0(indy(k), indx(k)) = npvsol0(indy(k), indx(k)) + La(k,4);
%     fnpsol0(findy(k), findx(k)) = fnpsol0(findy(k), findx(k)) + 1;
%     fnpusol0(findy(k), findx(k)) = fnpusol0(findy(k), findx(k)) + La(k,3);
%     fnpvsol0(findy(k), findx(k)) = fnpvsol0(findy(k), findx(k)) + La(k,4);
% end
% [npsol0, npusol0, npvsol0] = ProcEuData(ny, nx, npsol0, npusol0, npvsol0);
% [fnpsol0, fnpusol0, fnpvsol0] = ProcEuData(fny, fnx, fnpsol0, fnpusol0, fnpvsol0);

file_name = sprintf('./data/time%05d.mat', Nt);
load(file_name)

indx = ceil( ( La(:,1) - domain(1) )/hx );
indy = ceil( ( La(:,2)  - domain(3) )/hy );
findx = ceil( ( La(:,1) - domain(1) )/fhx );
findy = ceil( ( La(:,2)  - domain(3) )/fhy );

for k=1:np
    npsol(indy(k), indx(k)) = npsol(indy(k), indx(k)) + 1;
    npusol(indy(k), indx(k)) = npusol(indy(k), indx(k)) + La(k,3);
    npvsol(indy(k), indx(k)) = npvsol(indy(k), indx(k)) + La(k,4);
    fnpsol(findy(k), findx(k)) = fnpsol(findy(k), findx(k)) + 1;
    fnpusol(findy(k), findx(k)) = fnpusol(findy(k), findx(k)) + La(k,3);
    fnpvsol(findy(k), findx(k)) = fnpvsol(findy(k), findx(k)) + La(k,4);
end

[npsol, npusol, npvsol] = ProcEuData(ny, nx, npsol, npusol, npvsol);
[fnpsol, fnpusol, fnpvsol] = ProcEuData(fny, fnx, fnpsol, fnpusol, fnpvsol);
npsol00 = npsol; npusol00 = npusol; npvsol00 = npvsol; 
fnpsol00 = fnpsol; fnpusol00 = fnpusol; fnpvsol00 = fnpvsol; 

[xx,yy] = meshgrid(linspace(-pi,pi,nx), linspace(-pi,pi,nx));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)];

for j=2:ns
    npsol0 = npsol; npusol0 = npusol; npvsol0 = npvsol; 

    uox = exp(1i * nx_vec * kk) * (u_hat(:,j-1) .* transpose(rk(1,:)));
    uoy = exp(1i * nx_vec * kk) * (u_hat(:,j-1) .* transpose(rk(2,:)));
    uox = real(reshape(uox, ny, nx));
    uoy = real(reshape(uoy, ny, nx));
    forcex = beta*(uox.*npsol0 - npusol0);
    forcey = beta*(uoy.*npsol0 - npvsol0);
    
    npsol = solveNumDensity(domain, nx, ny, dt, npusol0./npsol0, npvsol0./npsol0, 0.0*npsol0, npsol0, 0.0);
    npusol = solveNumDensity(domain, nx, ny, dt, npusol0./npsol0, npvsol0./npsol0, forcex, npusol0, 0.0);
    npvsol = solveNumDensity(domain, nx, ny, dt, npusol0./npsol0, npvsol0./npsol0, forcey, npvsol0, 0.0);
end

[xx,yy] = meshgrid(linspace(-pi,pi,fnx), linspace(-pi,pi,fnx));
nx_vec = [reshape(xx,[],1), reshape(yy,[],1)];

for j=2:ns
    fnpsol0 = fnpsol; fnpusol0 = fnpusol; fnpvsol0 = fnpvsol; 

    uox = exp(1i * nx_vec * kk) * (u_hat(:,j-1) .* transpose(rk(1,:)));
    uoy = exp(1i * nx_vec * kk) * (u_hat(:,j-1) .* transpose(rk(2,:)));
    uox = real(reshape(uox, fny, fnx));
    uoy = real(reshape(uoy, fny, fnx));
    forcex = beta*(uox.*fnpsol0 - fnpusol0);
    forcey = beta*(uoy.*fnpsol0 - fnpvsol0);

    fnpsol = solveNumDensity(domain, fnx, fny, dt, fnpusol0./fnpsol0, fnpvsol0./fnpsol0, 0.0*fnpsol0, fnpsol0, 0.0);
    fnpusol = solveNumDensity(domain, fnx, fny, dt, fnpusol0./fnpsol0, fnpvsol0./fnpsol0, forcex, fnpusol0, 0.0);
    fnpvsol = solveNumDensity(domain, fnx, fny, dt, fnpusol0./fnpsol0, fnpvsol0./fnpsol0, forcey, fnpvsol0, 0.0);
end



dfu = fnpusol - fnpusol0;
dfv = fnpvsol - fnpvsol0;
du = npusol - npusol0;
dv = npvsol - npvsol0;
c2f = ones(n2); %/n2^2; % n2^2 here is to ensure <pho v> in the same scale when in fine mesh
duu = kron(du,c2f); dvv = kron(dv,c2f);

sigx = std(reshape(dfu - duu, [], 1));
sigy = std(reshape(dfv - dvv, [], 1));

sig = [sigx sigy];

% figure; 
% subplot(3,4,1); imagesc(fnpusol00); colorbar
% subplot(3,4,2); imagesc(fnpvsol00); colorbar
% subplot(3,4,3); imagesc(npusol00); colorbar
% subplot(3,4,4); imagesc(npvsol00); colorbar
% 
% % subplot(3,4,1); imagesc(fnpusol0); colorbar
% % subplot(3,4,2); imagesc(fnpvsol0); colorbar
% % subplot(3,4,3); imagesc(npusol0); colorbar
% % subplot(3,4,4); imagesc(npvsol0); colorbar
% 
% subplot(3,4,5); imagesc(fnpusol); colorbar
% subplot(3,4,6); imagesc(fnpvsol); colorbar
% subplot(3,4,7); imagesc(npusol); colorbar
% subplot(3,4,8); imagesc(npvsol); colorbar
% 
% subplot(3,4,9); imagesc(dfu); colorbar
% subplot(3,4,10); imagesc(dfv); colorbar
% subplot(3,4,11); imagesc(duu); colorbar
% subplot(3,4,12); imagesc(dvv); colorbar

end