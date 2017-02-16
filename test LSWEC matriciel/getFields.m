function [u,v,h]=getFields(lat,lon,lev,time,waveFlag)

        % number of grid points
        nj = length(lat);
        ni = length(lon);
        nl = length(lev);
        nt = length(time);

        % preallocating
        u = zeros(nj,ni,nl,nt);
        v = zeros(nj,ni,nl,nt);
        h = zeros(nj,ni,nl,nt);

        % phase speed
        C  = getPhaseSpeed(waveFlag);

        % latitude-dependent amplitudes
        [uTilde,vTilde,hTilde] = getAmplitudes(lat,waveFlag);

        % adding time and longitude dependence
        for t=1:nt
            for l=1:nl
                u(:,:,nl,t) = real(uTilde*exp(1i*k*(lon-C*t(t))));
                v(:,:,nl,t) = real(vTilde*exp(1i*k*(lon-C*t(t))));
                h(:,:,nl,t) = real(hTilde*exp(1i*k*(lon-C*t(t))));    
            end
        end

    end