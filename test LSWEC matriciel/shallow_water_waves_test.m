function [u,v,h,Ps,T]=shallow_water_waves_test(lat,lon,lev,time,waveFlag)
% SHALLOW_WATER_WAVES_TEST provides the analytic fields for the shallow
% water test case
% 
% This code is part of the online supporting material for: Shamir O, Paldor
% N. A quantitative test case for global-scale dynamical cores based on
% analytic wave solutions of the Shallow Water Equations. 
% Submitted to: Q J ROY METEOR SOC, Feb 2016.
% 
% The following functions calculate the analytic fields for the proposed
% test case on arbitrary lat x lon grids (See common arguments list below).
% 
% C = getPhaseSpeed(waveFlag):
%   Returns the analytic phase speed for the desired wave.
% 
% [u, v, h] = getFields(lat,lon,lev,time,waveFlag):
%   Returns the analytic velocity and free-surface height anomaly fields.
% 
% Ps = getSurfacePressure(h):
%   Transforms the free-surface height field into surface pressure.
% 
% T = getInitialTemperature(lat,lon,lev,time):
%   Returns the initial temperature field used in the text.
% 
% Arguments:
%     lat       - *column* vector of desired latitudes (radians) for output
%     lon       - *row* vector of desired longitudes (radians) for output
%     lev       - vector of desired vertical levels for output 
%     time      - vector of desired times (sec) for output
%     waveFlag  - WIG=-1, Rossby=0, EIG=1 
%     C         - phase speed (rad/sec) 
%     u         - zonal velocity field (m/sec) (lat,lon,lev,time) 
%     v         - meridional velocity field (m/sec) (lat,lon,lev,time)
%     h         - free-surface height anomaly field (m) (lat,lon,time)
%     Ps        - surface pressure (Pa) (lat,lon,time)
%     T         - initial temperature field (Kelvin) (lat,lon,time)    
% 
% Notes:
%   1)  Since the analytic velocity fields of the proposed test case are
%       vertically homogeneous, the choice of vertical levels in getFields
%       and getInitialTemperature is arbitrary. Only the length of the 1D
%       lev array is actually used.
% 
%   2)  In accordance with the proposed test's parameters, this code is
%       only valid for gH=5e4 m^2/s^2 and (n,k)=(5,10). Attempting to use
%       this code as is with different values would yield erroneous
%       analytic solutions.
%      
      
    % parameters
    omega = 7.29212e-5;          % Earth's angular frequency (rad/sec)
    g     = 10.0;                % Earth's gravitational acceleration (m/sec^2)
    a     = 6371220.0;           % Earth's mean radius (m)
    H0    = 5.0e3;               % Layer's mean depth (m)
    P0    = 1.0e5;               % Reference pressure (Pa)
    Rd    = 287.0;               % Gas const for dry air (J/kg/K)
    n     = 5;                   % chosen mode number
    k     = 10;                  % chosen wave number
    sigma = 0.5+(0.25+k^2)^0.5;  % see Eq.10 in text
    amp   = 1.e-8;               % arbitrary amplitude for linear waves
    
    [u,v,h] = getFields(lat,lon,lev,time,waveFlag);
    
    Ps = getSurfacePressure(h);
    
    T = getInitialTemperature(lat,lon,lev,time);

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

    function Ps=getSurfacePressure(h)
        Ps = P0*(1+h/H0);
    end

    function T=getInitialTemperature(lat,lon,lev,time)
	
        % number of grid points
        nj = length(lat);
        ni = length(lon);
        nl = length(lev);
        nt = length(time);

        T  = ones(nj,ni,nl,nt)*g*H0/Rd;
    end

end