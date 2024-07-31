 
sigma_x = 0.01;beta = 1; L = 20;D = 1;sigma_k = 0.05;
func = @(r) r-(-beta + sqrt(beta^2+2/sigma_x^2/(D+beta+1/sigma_x^2*r)*L*(-D + sqrt(D^2 + 1/sigma_x^2*(1/(D+beta+r/sigma_x^2))^2*L*sigma_k^2) )/((1/(D+beta+r/sigma_x^2))^2*L/sigma_x^2)  ))*sigma_x^2;
r0 = 1; % starting point
roots = fzero(func,r0) ;
alpha = 1/(D+beta+1/sigma_x^2*roots);
r2_rd = (-D + sqrt(D^2+1/sigma_x^2*alpha^2*L*sigma_k^2))*sigma_x^2/alpha^2/L;
% disp('Theoretic value of the diagonal entry of R22')
disp(r2_rd)
 

  