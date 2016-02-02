function [ g ] = atsf(lambda, teta,time,c,space,f,filtre,delta)
%% ************************************************************************
% author : Brachet Matthieu (IECL-Metz)
% birth date : janv-30-2016
% last version : 01-feb-2016
% *************************************************************************
% this function build the amplification number obtain with Runge Kutta method is time,
% space discretisation and an other scheme for space discretisation and a space
% filter on the advection equation.
%
% *** inputs
% lambda : cfl number;
% teta : wave number;
% time : type of time method ('rk1', 'rk2', 'rk4', 'dirk12','dirk23' or
% 'dirk34')
% c : order for the space scheme (2 (only for explicit),4 or 6).
% space : type of discretisation (explicit or implicit).
% f : order of the space filter(redonnet, long or shapiro);
% filtre : type of space filter (filter = 0 : explicit, filter = 1,
% implicit).
%
% *** outputs
% g : amplification number.
%
% *************************************************************************

%% *** calcul du filtre ***************************************************
if strcmp(filtre,'redonnet') == 1
if f==0
    ftr=1;
elseif f==2
    ftr = 0.5+2*(1/4*cos(teta));
elseif f==4
    ftr = 10/16+2*(4/16*cos(teta)-1/16*cos(2*teta));
elseif f==6
    ftr = 44/64+2*(15/64*cos(teta)-6/64*cos(2*teta)+1/64*cos(3*teta));
elseif f == 8
    ftr = 186/256+2*(56/256*cos(teta)-28/256*cos(2*teta)+8/256*cos(3*teta)-1/256*cos(4*teta));
elseif f == 10
    ftr = 772/1024+2*(210/1024*cos(teta)-120/1024*cos(2*teta)+45/1024*cos(3*teta)-10/1024*cos(4*teta)+1/1024*cos(5*teta));
else
    error('Invalid input argument : f must be an interger in {0,2,4,6,8,10} for explicit filter.');
end

elseif strcmp(filtre,'long') == 1
if f==0
    ftr=1;
elseif f==2
    ftr = 1./(1+delta.*tan(teta/2).^2);
else
    error('Invalid input argument : f must be an interger in {0,2} for implicit filter.');
end

elseif strcmp(filtre,'shapiro') == 1
    ftr=1-sin(teta/2).^f;
end

%% *** space discretisation ***********************************************
if strcmp(space,'explicit') == 1
    
if c == 2
    s=sin(teta);
elseif c == 4
    s=4/3*sin(teta)-1/3*sin(2*teta);
elseif c == 6
    s=3/2*sin(teta)-3/5*sin(2*teta)+1/10*sin(3*teta);
else
    error('Invalid input argument : c must be an integer in {2,4,6} for explicit space scheme.');
end

elseif strcmp(space,'implicit') == 1

if c == 4
    s = sin(teta)./(2/3+2*(1/6*cos(teta)));
elseif c == 6
    s = (4/3*sin(teta)-1/3*sin(2*teta))./(4/5+2*(2/15*cos(teta)-1/30*cos(2*teta)));
else
    error('Invalid input argument : c must be an integer in {4,6} for implicit space scheme.');
end

end

%% *** time discretisation (only with rk4) ********************************
% be careful : here i is the complex number such that 1i*1i = -1 !
if strcmp(time,'rk4') == 1
    g = 1-i*lambda.*s-1/2*lambda.^2.*s.^2+1/6*i*lambda.^3.*s.^3+1/24*lambda.^4.*s.^4;
    g=g.*ftr;
elseif strcmp(time,'rk1') == 1
    g=1-i*lambda.*s;
    g=g.*ftr;
elseif strcmp(time,'rk2') == 1
    g = 1-i*lambda.*s-1/2*lambda.^2.*s.^2;
    g=g.*ftr;
elseif strcmp(time,'dirk12') == 1
    theta = -i*lambda*s;
    g = (1+theta/2)./(1-theta/2);
    g = g.*ftr;
elseif strcmp(time,'dirk23') == 1
    theta = -i.*lambda.*s;
    nom = -2.*(18-9*theta-9*sqrt(3)*theta+3*theta.^3+2*sqrt(3).*theta.^3);
    denom = (-6 + 3.*theta + theta.*sqrt(3)).*(6-6*theta-2*sqrt(3).*theta+2*theta.^2+sqrt(3).*theta.^2);
    g=nom./denom;
    g = g.*ftr;
elseif strcmp(time,'dirk34') == 1
    theta=-i.*lambda.*s;
    alfa=2*cos(pi/18)/sqrt(3);
    num = -24*alfa+12*theta*alfa+36*theta*alfa.^2 +6.*theta.^2*alfa -18*theta.^2.*alfa^3 - 7*theta.^3.*alfa-9*theta.^3*alfa^3+3*theta.^3*alfa^4-2*theta.^3;
    denom = 3.*alfa.*(-2+theta+theta.*alfa).^3;
    g=num./denom;
    g=g.*ftr;
else
    error('Invalid input argument : time scheme (in time) is invalid');
end
    
    
end

