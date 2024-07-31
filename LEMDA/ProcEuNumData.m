function [sol] = ProcEuNumData(ny, nx, npsol)
mx = max(npsol(:)); mx = 1.0/mx/nx/ny;

if (npsol(1,1) < 1)
    npsol(1,1) = mx;
end

if (npsol(1,nx) < 1)
    npsol(1,nx) = mx;
end

if (npsol(ny,1) < 1)
    npsol(ny,1) = mx;
end

if (npsol(ny,nx) < 1)
    npsol(ny,nx) = mx;
end

for jy = 2:ny-1
    jx = 1;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
    end

    jx = nx;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
    end
end

for jx = 2:nx-1
    jy = 1;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
    end

    jy = ny;
    if (npsol(jy,jx) < 1)
        npsol(jy,jx) = mx;
    end
end

for jy = 2:ny-1
    for jx = 2:nx-1
        if (npsol(jy,jx) < 1)
            npsol(jy,jx) = mx;
        end
    end
end

sol = npsol;
end