function usol = solveNumDensity(domain, nx, ny, dt, vx, vy, force, sol, sigma_n)

rhs = zeros(size(sol));
% hx = ( domain(2) - domain(1) )./(nx-1);
% hy = ( domain(4) - domain(3) )./(ny-1);
hx = ( domain(2) - domain(1) )./nx;
hy = ( domain(4) - domain(3) )./ny;

%ssol = sol;
for j=1:ny
    for k=1:nx
        tem = 0.0;
        % finite volume fluxes; essentially the same as the center FDM here
%         tem = tem + 0.5*(vx(j,k)*sol(j,k) + vx(j,k+1)*sol(j,k+1))/hx; % right boundary of a control volume/cell
%         tem = tem - 0.5*(vx(j,k)*sol(j,k) + vx(j,k-1)*sol(j,k-1))/hx; % left boundary
%         tem = tem - 0.5*(vy(j,k)*sol(j,k) + vy(j-1,k)*sol(j-1,k))/hy; % bottom boundary
%         tem = tem + 0.5*(vy(j,k)*sol(j,k) + vy(j+1,k)*sol(j+1,k))/hy; % top boundary

%         % upwind finite difference scheme
%         if (vx(j,k)>0 && k>1)
%             tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
%         elseif (vx(j,k)>0 && k==1)
%             tem = tem + (vx(j,k)*sol(j,k) - vx(j,nx)*sol(j,nx)) / hx;
%         elseif (vx(j,k)<0 && k<nx)
%             tem = tem - (vx(j,k)*sol(j,k) - vx(j,k+1)*sol(j,k+1)) / hx;
%         elseif (vx(j,k)<0 && k==nx)
%             tem = tem - (vx(j,k)*sol(j,k) - vx(j,1)*sol(j,1)) / hx;
%         end
% 
%         if (vy(j,k)>0 && j>1)
%             tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         elseif (vy(j,k)>0 && j==1)
%             tem = tem + (vy(j,k)*sol(j,k) - vy(ny,k)*sol(ny,k)) / hy;
%         elseif (vy(j,k)<0 && j<ny)
%             tem = tem - (vy(j,k)*sol(j,k) - vy(j+1,k)*sol(j+1,k)) / hy;
%         elseif (vy(j,k)<0 && j==ny)
%             tem = tem - (vy(j,k)*sol(j,k) - vy(1,k)*sol(1,k)) / hy;
%         end
                
        % Lax-Wendroff scheme
        % tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx &
        %      - 0.5*dt*( 0.5*(vx(j,k+1) + vx(j,k))*(sol(j,k+1)-sol(j,k)) - 0.5*(vx(j,k-1) + vx(j,k))*(sol(j,k) - sol(j,k-1)) ) /hx/hx
        % tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy &
        %      - 0.5*dt*( 0.5*(vy(j+1,k) + vx(j,k))*(sol(j+1,k)-sol(j,k)) - 0.5*(vx(j-1,k) + vy(j,k))*(sol(j,k) - sol(j-1,k)) ) /hy/hy

        % Lax-Friedrichs scheme (with periodic boundary conditions)
        if (j==1 && k==1)
            %tem = tem + (vx(j,k+1)*sol(j,k+1) - vx(j,k)*sol(j,k)) / hx;
            %tem = tem + (vy(j+1,k)*sol(j+1,k) - vy(j,k)*sol(j,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j,k+1) + 2*sol(j,k) );
            
            tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,nx)*sol(j,nx)) / hx;
            tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(ny,k)*sol(ny,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(ny,k) + sol(j,k+1) + sol(j,nx) );
        elseif (j==1 && k==nx)
            %tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
            %tem = tem + (vy(j+1,k)*sol(j+1,k) - vy(j,k)*sol(j,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + 2*sol(j,k) + sol(j,k-1) );
            
            tem = tem + 0.5*(vx(j,1)*sol(j,1) - vx(j,k-1)*sol(j,k-1)) / hx;
            tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(ny,k)*sol(ny,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(ny,k) + sol(j,1) + sol(j,k-1) );
        elseif (j==ny && k==nx)
            %tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
            %tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j-1,k) + 2*sol(j,k) + sol(j,k-1) );
            
            tem = tem + 0.5*(vx(j,1)*sol(j,1) - vx(j,k-1)*sol(j,k-1)) / hx;
            tem = tem + 0.5*(vy(1,k)*sol(1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(1,k) + sol(j-1,k) + sol(j,1) + sol(j,k-1) );
        elseif (j==ny && k==1)
            %tem = tem + (vx(j,k+1)*sol(j,k+1) - vx(j,k)*sol(j,k)) / hx;
            %tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j-1,k) + sol(j,k+1) + 2*sol(j,k) );
            
            tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,nx)*sol(j,nx)) / hx;
            tem = tem + 0.5*(vy(1,k)*sol(1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(1,k) + sol(j-1,k) + sol(j,k+1) + sol(j,nx) );
        elseif (j==1)
            %tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
            %tem = tem + (vy(j+1,k)*sol(j+1,k) - vy(j,k)*sol(j,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j,k) + sol(j,k+1) + sol(j,k-1) );
            
            tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
            tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(ny,k)*sol(ny,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(ny,k) + sol(j,k+1) + sol(j,k-1) );
        elseif (j==ny)
            %tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
            %tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j,k) + sol(j-1,k) + sol(j,k+1) + sol(j,k-1) );
            
            tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
            tem = tem + 0.5*(vy(1,k)*sol(1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(1,k) + sol(j-1,k) + sol(j,k+1) + sol(j,k-1) );
        elseif (k==nx)
            %tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
            %tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j-1,k) + sol(j,k) + sol(j,k-1) );
            
            tem = tem + 0.5*(vx(j,1)*sol(j,1) - vx(j,k-1)*sol(j,k-1)) / hx;
            tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j-1,k) + sol(j,1) + sol(j,k-1) );
        elseif (k==1)
            %tem = tem + (vx(j,k+1)*sol(j,k+1) - vx(j,k)*sol(j,k)) / hx;
            %tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j-1,k) + sol(j,k+1) + sol(j,k) );
            
            tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,nx)*sol(j,nx)) / hx;
            tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j-1,k) + sol(j,k+1) + sol(j,nx) );
        else
            tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
            tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
            %ssol(j,k) = 0.25*( sol(j+1,k) + sol(j-1,k) + sol(j,k+1) + sol(j,k-1) );
        end

        rhs(j,k) = -tem + force(j,k);
    end
end

% Two-step Adams–Bashforth
%usol = ssol + 1.5*dt*rhs - 0.5*rhs0*dt

%usol = ssol + dt*rhs + randn(size(sol)) * sigma_n * sqrt(dt); % Lax-Friedrichs scheme; too disspasive w
usol = sol + dt*rhs + randn(size(sol)) * sigma_n * sqrt(dt); % center-difference
end
