%if(indx>1 && indx<nx && indy>1 && indy <ny)
if (indx>1)
    indk = (indy-1)*nx + indx-1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (indy-1)*nx + nx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)-a y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end

if (indx<nx)
    indk = (indy-1)*nx + indx+1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (indy-1)*nx + 1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)+a y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end

if (indy>1)
    indk = (indy-2)*nx + indx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (ny-1)*nx + indx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki)-b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end

if (indy<ny)
    indk = (indy-0)*nx + indx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (0)*nx + indx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki)+b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end

if (indx==1 && indy==1) % left bottom corner
    indk = ny*nx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)-a y(ki)-b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx==1 && indy>1)
    indk = (indy-1)*nx; % (indy-2)*nx + nx
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)-a y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx>1 && indy==1)
    indk = (ny-1)*nx + indx-1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki)-b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (indy-2)*nx + indx-1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end

if (indx==nx && indy==1) % right bottom corner
    indk = (ny-1)*nx + 1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)+a y(ki)-b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx==nx && indy>1)
    indk = (indy-2)*nx + 1; 
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)+a y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx<nx && indy==1)
    indk = (ny-1)*nx + indx+1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki)-b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (indy-2)*nx + indx+1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end


if (indx==1 && indy==ny) % top left corner
    indk = nx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)-a y(ki)+b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx==1 && indy<ny)
    indk = (indy)*nx + nx;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)-a y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx>1 && indy==ny)
    indk = (0)*nx + indx-1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki)-b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (indy-0)*nx + indx-1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end

if (indx==nx && indy==ny) % top right corner
    indk = 1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)+a y(ki)+b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx==nx && indy<ny)
    indk = (indy)*nx + 1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki)+a y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
elseif (indx<nx && indy==ny)
    indk = (0)*nx + indx+1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki)+b u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
else
    indk = (indy-0)*nx + indx+1;
    for k=1:cell4p(indk, 1)
        ki = cell4p(indk, k+1);
        floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
        [CFtem, ~] = getCF(floej, floek);
        CFx(j) = CFx(j) + CFtem(1);
        CFy(j) = CFy(j) + CFtem(2);
    end
end


%         if(indx>1 && indx<nx && indy>1 && indy <ny)
%             for jj=0:2
%                 indk = (indy-jj)*nx + indx-1;
%                 for k=1:cell4p(indk, 1)
%                     ki = cell4p(indk, k+1);
%                     floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
%                     [CFtem, ~] = getCF(floej, floek);
%                     CFx(j) = CFx(j) + CFtem(1);
%                     CFy(j) = CFy(j) + CFtem(2);
%                 end
%
%                 indk = (indy-jj)*nx + indx+1;
%                 for k=1:cell4p(indk, 1)
%                     ki = cell4p(indk, k+1);
%                     floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
%                     [CFtem, ~] = getCF(floej, floek);
%                     CFx(j) = CFx(j) + CFtem(1);
%                     CFy(j) = CFy(j) + CFtem(2);
%                 end
%             end
%
%             indk = (indy-2)*nx + indx;
%             for k=1:cell4p(indk, 1)
%                 ki = cell4p(indk, k+1);
%                 floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
%                 [CFtem, ~] = getCF(floej, floek);
%                 CFx(j) = CFx(j) + CFtem(1);
%                 CFy(j) = CFy(j) + CFtem(2);
%             end
%
%             indk = (indy-0)*nx + indx;
%             for k=1:cell4p(indk, 1)
%                 ki = cell4p(indk, k+1);
%                 floek = [r(ki) x(ki) y(ki) u(ki) v(ki) 0.0];
%                 [CFtem, ~] = getCF(floej, floek);
%                 CFx(j) = CFx(j) + CFtem(1);
%                 CFy(j) = CFy(j) + CFtem(2);
%             end
%         end
