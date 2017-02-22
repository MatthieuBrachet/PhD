function [u,v,h]=getAmplitudes(lat,waveFlag)

        cp = cos(lat); 
        sp = sin(lat);
        tp = tan(lat);

        C  = getPhaseSpeed(waveFlag);

        [psi,dpsi] = getPsi(lat);

        switch waveFlag
            case 'Rossby'
                Kp = (g*H0+(C*a*cp).^2)./(C*cp); 
                Km = (g*H0-(C*a*cp).^2)./(C*cp); 
                V  = sqrt(2*omega*abs(Km))./cp.*psi;
                h  = sqrt(2*omega*a^2*H0^2*abs(Km))./Km.*...
                     (dpsi+tp.*(0.5*Kp./Km-2*omega/C).*psi);
                u  = 2*omega*sp/C.*V + g./(a*C*cp).*h;
            case {'EIG','WIG'}
                W1 = ((k*C)^2-(2*omega*sp).^2)./(C*cp); 
                W2 = ((k*C)^2-(2*omega*sp).^2)./(2*omega*cp.^2); 
                V  = sqrt(2*omega*g*H0*abs(W1))./(W1.*cp).*...
                     (dpsi+tp.*(0.5-2*omega./W2+2*omega/C).*psi);
                h  = sqrt(2*omega*a^2*H0/g*abs(W1)).*psi;
                u  = 2*omega*sp/C.*V + g./(a*C*cp).*h;
        end

        v      = -1i*k*V;

    end