function [cff, flag] = getCF(floei, floej)
    flag = 0;
    cff = zeros(1, 4);

    ri = floei(1);
    xi = floei(2:3);
    vi = floei(4:6);

    rj = floej(1);
    xj = floej(2:3);
    vj = floej(4:6);

    if ri < rj
        rmax = rj;
        rmin = ri;
    else
        rmax = ri;
        rmin = rj;
    end

    Eij = 2.0e4;
    Gij = 2.0e4;
    muij = 0.2;

    nx = xi(1) - xj(1);
    ny = xi(2) - xj(2);
    distij = sqrt(nx^2 + ny^2);
    delta = distij - ri - rj;

    if delta < 0.0
        flag = 1;
        nn(1) = nx / distij;
        nn(2) = ny / distij;

        tem = 4.0 * distij^2 * rmax^2 - (distij^2 - rmin^2 + rmax^2)^2;
        cl = sqrt(abs(tem)) / distij;

        ncf(1) = -cl * Eij * delta * nn(1);
        ncf(2) = -cl * Eij * delta * nn(2);

        vd(1) = vj(1) - vi(1);
        vd(2) = vj(2) - vi(2);
        avd(1) = -rj * nn(2) * vj(3) - ri * nn(2) * vi(3);
        avd(2) = rj * nn(1) * vj(3) + ri * nn(1) * vi(3);
        vt = -(vd(1) + avd(1)) * nn(2) + (vd(2) + avd(2)) * nn(1);
        tcf(1) = -cl * Gij * vt * nn(2);
        tcf(2) = cl * Gij * vt * nn(1);

        a = sqrt(tcf(1)^2 + tcf(2)^2);
        b = muij * sqrt(ncf(1)^2 + ncf(2)^2);

        if a > b % Coulomb's Law of Friction 
            tcf = tcf * b / a;
        end

        cff(1) = ncf(1) + tcf(1);
        cff(2) = ncf(2) + tcf(2);
        cff(3) = ri * (nn(1) * tcf(2) - nn(2) * tcf(1));
        cff(4) = rj * (nn(1) * tcf(2) - nn(2) * tcf(1));

        if distij < 1.0e-7
            disp('Concentric floes detected!');
            cff = zeros(1, 4);
        end

        if any(isnan(cff))
            disp('Abnormal collision detected:');
            disp(['distij: ' num2str(distij) ' tem: ' num2str(tem)]);
            disp(['floei: ' num2str(floei)]);
            disp(['floej: ' num2str(floej)]);
            cff = zeros(1, 4);
        end
    end
end