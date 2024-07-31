function [sol, usol, vsol] = ProcEuData(ny, nx, npsol, npusol, npvsol)
mx = max(npsol(:)); mx = 1.0/mx/nx/ny;

if (npsol(1,1) < 1)
    npsol(1,1) = mx;
%     npsol(1,1) = (npsol(2,1) + npsol(1,2) + npsol(2,2))/3;
%     npusol(1,1) = (npusol(2,1) + npusol(1,2) + npusol(2,2))/3;
%     npvsol(1,1) = (npvsol(2,1) + npvsol(1,2) + npvsol(2,2))/3;
end

if (npsol(1,nx) < 1)
    npsol(1,nx) = mx;
%     npsol(1,nx) = (npsol(2,nx) + npsol(1,nx-1) + npsol(2,nx-1))/3;
%     npusol(1,nx) = (npusol(2,nx) + npusol(1,nx-1) + npusol(2,nx-1))/3;
%     npvsol(1,nx) = (npvsol(2,nx) + npvsol(1,nx-1) + npvsol(2,nx-1))/3;
end

if (npsol(ny,1) < 1)
    npsol(ny,1) = mx;
%     npsol(ny,1) = (npsol(ny-1,1) + npsol(ny,2) + npsol(ny-1,2))/3;
%     npusol(ny,1) = (npusol(ny-1,1) + npusol(ny,2) + npusol(ny-1,2))/3;
%     npvsol(ny,1) = (npvsol(ny-1,1) + npvsol(ny,2) + npvsol(ny-1,2))/3;
end

if (npsol(ny,nx) < 1)
    npsol(ny,nx) = mx;
%     npsol(ny,nx) = (npsol(ny-1,nx) + npsol(ny,nx-1) + npsol(ny-1,nx-1))/3;
%     npusol(ny,nx) = (npusol(ny-1,nx) + npusol(ny,nx-1) + npusol(ny-1,nx-1))/3;
%     npvsol(ny,nx) = (npvsol(ny-1,nx) + npvsol(ny,nx-1) + npvsol(ny-1,nx-1))/3;
end

for jy = 2:ny-1
    jx = 1;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
%         npsol(jy,jx) = ( npsol(jy-1,jx) + npsol(jy+1,jx) + npsol(jy,jx+1) )/3;
%         npusol(jy,jx) = ( npusol(jy-1,jx) + npusol(jy+1,jx) + npusol(jy,jx+1) )/3;
%         npvsol(jy,jx) = ( npvsol(jy-1,jx) + npvsol(jy+1,jx) + npvsol(jy,jx+1) )/3;
    end

    jx = nx;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
%         npsol(jy,jx) = ( npsol(jy-1,jx) + npsol(jy,jx-1) + npsol(jy+1,jx) )/3;
%         npusol(jy,jx) = ( npusol(jy-1,jx) + npusol(jy,jx-1) + npusol(jy+1,jx) )/3;
%         npvsol(jy,jx) = ( npvsol(jy-1,jx) + npvsol(jy,jx-1) + npvsol(jy+1,jx) )/3;
    end
end

for jx = 2:nx-1
    jy = 1;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
%         npsol(jy,jx) = ( npsol(jy,jx-1) + npsol(jy+1,jx) + npsol(jy,jx+1) )/3;
%         npusol(jy,jx) = ( npusol(jy,jx-1) + npusol(jy+1,jx) + npusol(jy,jx+1) )/3;
%         npvsol(jy,jx) = ( npvsol(jy,jx-1) + npvsol(jy+1,jx) + npvsol(jy,jx+1) )/3;
    end

    jy = ny;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
%         npsol(jy,jx) = ( npsol(jy-1,jx) + npsol(jy,jx-1) + npsol(jy,jx+1) )/3;
%         npusol(jy,jx) = ( npusol(jy-1,jx) + npusol(jy,jx-1) + npusol(jy,jx+1) )/3;
%         npvsol(jy,jx) = ( npvsol(jy-1,jx) + npvsol(jy,jx-1) + npvsol(jy,jx+1) )/3;
    end
end

for jy = 2:ny-1
    for jx = 2:nx-1
        if (npsol(jy,jx) < 1)
            npsol(jy,jx) = mx;
%             npsol(jy,jx) = ( npsol(jy-1,jx) + npsol(jy,jx-1) + npsol(jy+1,jx) + npsol(jy,jx+1) )/4;
%             npusol(jy,jx) = ( npusol(jy-1,jx) + npusol(jy,jx-1) + npusol(jy+1,jx) + npusol(jy,jx+1) )/4;
%             npvsol(jy,jx) = ( npvsol(jy-1,jx) + npvsol(jy,jx-1) + npvsol(jy+1,jx) + npvsol(jy,jx+1) )/4;
        end
    end
end

sol = npsol; usol = npusol; vsol = npvsol;
end