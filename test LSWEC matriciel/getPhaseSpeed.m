function C=getPhaseSpeed(waveFlag)

        En     = g*H0/a^2*(n+sigma)^2;
        delta0 = 3*k^2*En;
        delta4 = -54*g*H0*omega*k^4/a^2;

        Cj = zeros(1,3);

        for j=1:3
            deltaj = (0.5*(delta4+sqrt(delta4^2-4*delta0^3)))^(1/3)*exp(2*pi*1i*j/3);
            Cj(j)  = -1/(3*k^2)*deltaj*(1+delta0/deltaj^2);
        end

        Cj = real(Cj);

        switch waveFlag
            case 'Rossby'
                C = -min(abs(Cj));
            case 'WIG'      
                C = min(Cj);
            case 'EIG'
                C = max(Cj);
        end

    end   