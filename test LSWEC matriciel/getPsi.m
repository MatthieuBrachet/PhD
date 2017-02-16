function [psi,dpsi]=getPsi(lat)

        sp = sin(lat);
        cp = cos(lat);
        tp = tan(lat);

        a3 = sigma*(sigma+1)*(sigma+2);
        a4 = a3*(sigma+3);
        a5 = a4*(sigma+4);

        G5  = (4*a5*sp.^5 - 20*a4*sp.^3 + 15*a3*sp )/15;
        G5p = (4*a5*sp.^4 - 12*a4*sp.^2 + 3*a3 ).*cp/3;

        psi  = amp*cp.^sigma.*G5;
        dpsi = amp*cp.^(sigma).*(-tp.*G5 + G5p);  

    end   