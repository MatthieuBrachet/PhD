function T=getInitialTemperature(lat,lon,lev,time)
	
        % number of grid points
        nj = length(lat);
        ni = length(lon);
        nl = length(lev);
        nt = length(time);

        T  = ones(nj,ni,nl,nt)*g*H0/Rd;
    end