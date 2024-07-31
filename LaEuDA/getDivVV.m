function dvv = getDivVV(domain, nx, ny, vx, vy, sol)

dvv = 0.0*sol;
% hx = ( domain(2) - domain(1) )./(nx-1);
% hy = ( domain(4) - domain(3) )./(ny-1);
hx = ( domain(2) - domain(1) )./nx;
hy = ( domain(4) - domain(3) )./ny;

for j=1:ny
    for k=1:nx
        tem = 0.0;

%         % Lax-Friedrichs scheme (with periodic boundary conditions)
%         if (j==1 && k==1)
%             %tem = tem + (vx(j,k+1)*sol(j,k+1) - vx(j,k)*sol(j,k)) / hx;
%             %tem = tem + (vy(j+1,k)*sol(j+1,k) - vy(j,k)*sol(j,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,nx)*sol(j,nx)) / hx;
%             tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(ny,k)*sol(ny,k)) / hy;
%         elseif (j==1 && k==nx)
%             %tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
%             %tem = tem + (vy(j+1,k)*sol(j+1,k) - vy(j,k)*sol(j,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,1)*sol(j,1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(ny,k)*sol(ny,k)) / hy;
%         elseif (j==ny && k==nx)
%             %tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
%             %tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,1)*sol(j,1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             tem = tem + 0.5*(vy(1,k)*sol(1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         elseif (j==ny && k==1)
%             %tem = tem + (vx(j,k+1)*sol(j,k+1) - vx(j,k)*sol(j,k)) / hx;
%             %tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,nx)*sol(j,nx)) / hx;
%             tem = tem + 0.5*(vy(1,k)*sol(1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         elseif (j==1)
%             %tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             %tem = tem + (vy(j+1,k)*sol(j+1,k) - vy(j,k)*sol(j,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(ny,k)*sol(ny,k)) / hy;
%         elseif (j==ny)
%             %tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             %tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             tem = tem + 0.5*(vy(1,k)*sol(1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         elseif (k==nx)
%             %tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
%             %tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,1)*sol(j,1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         elseif (k==1)
%             %tem = tem + (vx(j,k+1)*sol(j,k+1) - vx(j,k)*sol(j,k)) / hx;
%             %tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%             
%             tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,nx)*sol(j,nx)) / hx;
%             tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         else
%             tem = tem + 0.5*(vx(j,k+1)*sol(j,k+1) - vx(j,k-1)*sol(j,k-1)) / hx;
%             tem = tem + 0.5*(vy(j+1,k)*sol(j+1,k) - vy(j-1,k)*sol(j-1,k)) / hy;
%         end

%         % upwind finite difference scheme
        if (vx(j,k)>0 && k>1)
            tem = tem + (vx(j,k)*sol(j,k) - vx(j,k-1)*sol(j,k-1)) / hx;
        elseif (vx(j,k)>0 && k==1)
            tem = tem + (vx(j,k)*sol(j,k) - vx(j,nx)*sol(j,nx)) / hx;
        elseif (vx(j,k)<0 && k<nx)
            tem = tem - (vx(j,k)*sol(j,k) - vx(j,k+1)*sol(j,k+1)) / hx;
        elseif (vx(j,k)<0 && k==nx)
            tem = tem - (vx(j,k)*sol(j,k) - vx(j,1)*sol(j,1)) / hx;
        end

        if (vy(j,k)>0 && j>1)
            tem = tem + (vy(j,k)*sol(j,k) - vy(j-1,k)*sol(j-1,k)) / hy;
        elseif (vy(j,k)>0 && j==1)
            tem = tem + (vy(j,k)*sol(j,k) - vy(ny,k)*sol(ny,k)) / hy;
        elseif (vy(j,k)<0 && j<ny)
            tem = tem - (vy(j,k)*sol(j,k) - vy(j+1,k)*sol(j+1,k)) / hy;
        elseif (vy(j,k)<0 && j==ny)
            tem = tem - (vy(j,k)*sol(j,k) - vy(1,k)*sol(1,k)) / hy;
        end

        dvv(j,k) = tem;
    end
end

end
