% data processing to fix when one cell doe
% 
for j=1:N
    if (npsol(1,1,j) < 1)
        npsol(1,1,j) = (npsol(2,1,j) + npsol(1,2,j) + npsol(2,2,j))/3;
        npusol(1,1,j) = (npusol(2,1,j) + npusol(1,2,j) + npusol(2,2,j))/3;
        npvsol(1,1,j) = (npvsol(2,1,j) + npvsol(1,2,j) + npvsol(2,2,j))/3; 

%         npsol(1,1,j) = (npsol(2,1,j) + npsol(1,2,j) + npsol(2,2,j) + npsol(1,3,j) + npsol(3,1,j) )/5;
%         npusol(1,1,j) = (npusol(2,1,j) + npusol(1,2,j) + npusol(2,2,j) + npusol(1,3,j) + npusol(3,1,j) )/5;
%         npvsol(1,1,j) = (npvsol(2,1,j) + npvsol(1,2,j) + npvsol(2,2,j) + npvsol(1,3,j) + npvsol(3,1,j) )/5; 
    end

    if (npsol(1,nx,j) < 1)
        npsol(1,nx,j) = (npsol(2,nx,j) + npsol(1,nx-1,j) + npsol(2,nx-1,j))/3;
        npusol(1,nx,j) = (npusol(2,nx,j) + npusol(1,nx-1,j) + npusol(2,nx-1,j))/3;
        npvsol(1,nx,j) = (npvsol(2,nx,j) + npvsol(1,nx-1,j) + npvsol(2,nx-1,j))/3;
    end

    if (npsol(ny,1,j) < 1)
        npsol(ny,1,j) = (npsol(ny-1,1,j) + npsol(ny,2,j) + npsol(ny-1,2,j))/3;
        npusol(ny,1,j) = (npusol(ny-1,1,j) + npusol(ny,2,j) + npusol(ny-1,2,j))/3;
        npvsol(ny,1,j) = (npvsol(ny-1,1,j) + npvsol(ny,2,j) + npvsol(ny-1,2,j))/3;     
    end

    if (npsol(ny,nx,j) < 1)
        npsol(ny,nx,j) = (npsol(ny-1,nx,j) + npsol(ny,nx-1,j) + npsol(ny-1,nx-1,j))/3;
        npusol(ny,nx,j) = (npusol(ny-1,nx,j) + npusol(ny,nx-1,j) + npusol(ny-1,nx-1,j))/3;
        npvsol(ny,nx,j) = (npvsol(ny-1,nx,j) + npvsol(ny,nx-1,j) + npvsol(ny-1,nx-1,j))/3;
    end

    for jy = 2:ny-1
        jx = 1;
        if (npsol(jy,jx,j) < 1)
            npsol(jy,jx,j) = ( npsol(jy-1,jx,j) + npsol(jy+1,jx,j) + npsol(jy,jx+1,j) )/3;
            npusol(jy,jx,j) = ( npusol(jy-1,jx,j) + npusol(jy+1,jx,j) + npusol(jy,jx+1,j) )/3;
            npvsol(jy,jx,j) = ( npvsol(jy-1,jx,j) + npvsol(jy+1,jx,j) + npvsol(jy,jx+1,j) )/3;
        end

        jx = nx;
        if (npsol(jy,jx,j) < 1)
            npsol(jy,jx,j) = ( npsol(jy-1,jx,j) + npsol(jy,jx-1,j) + npsol(jy+1,jx,j) )/3;
            npusol(jy,jx,j) = ( npusol(jy-1,jx,j) + npusol(jy,jx-1,j) + npusol(jy+1,jx,j) )/3;
            npvsol(jy,jx,j) = ( npvsol(jy-1,jx,j) + npvsol(jy,jx-1,j) + npvsol(jy+1,jx,j) )/3;
        end
    end

    for jx = 2:nx-1
        jy = 1;
        if (npsol(jy,jx,j) < 1)
            npsol(jy,jx,j) = ( npsol(jy,jx-1,j) + npsol(jy+1,jx,j) + npsol(jy,jx+1,j) )/3;
            npusol(jy,jx,j) = ( npusol(jy,jx-1,j) + npusol(jy+1,jx,j) + npusol(jy,jx+1,j) )/3;
            npvsol(jy,jx,j) = ( npvsol(jy,jx-1,j) + npvsol(jy+1,jx,j) + npvsol(jy,jx+1,j) )/3;
        end

        jy = ny;
        if (npsol(jy,jx,j) < 1)
            npsol(jy,jx,j) = ( npsol(jy-1,jx,j) + npsol(jy,jx-1,j) + npsol(jy,jx+1,j) )/3;
            npusol(jy,jx,j) = ( npusol(jy-1,jx,j) + npusol(jy,jx-1,j) + npusol(jy,jx+1,j) )/3;
            npvsol(jy,jx,j) = ( npvsol(jy-1,jx,j) + npvsol(jy,jx-1,j) + npvsol(jy,jx+1,j) )/3;
        end
    end

    for jy = 2:ny-1
        for jx = 2:nx-1
            if (npsol(jy,jx,j) < 1)
                npsol(jy,jx,j) = ( npsol(jy-1,jx,j) + npsol(jy,jx-1,j) + npsol(jy+1,jx,j) + npsol(jy,jx+1,j) )/4;
                npusol(jy,jx,j) = ( npusol(jy-1,jx,j) + npusol(jy,jx-1,j) + npusol(jy+1,jx,j) + npusol(jy,jx+1,j) )/4;
                npvsol(jy,jx,j) = ( npvsol(jy-1,jx,j) + npvsol(jy,jx-1,j) + npvsol(jy+1,jx,j) + npvsol(jy,jx+1,j) )/4;
            end
        end
    end
end
    