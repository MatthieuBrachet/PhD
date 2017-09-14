function [ g ] = atsf(lambda, teta,time,c,space,f,filtre,delta)
%% ************************************************************************
% author : Brachet Matthieu (IECL-Metz)
% birth date : janv-30-2016
% last version : 09-feb-2016
% *************************************************************************
% this function build the amplification number obtain with Runge Kutta method is time,
% space discretisation and an other scheme for space discretisation and a space
% filter on the advection equation.
%
% *** inputs
% lambda : cfl number;
% teta : wave number;
% time : type of time method ('rk1', 'rk2', 'rk4', 'dirk12a','dirk23a' or
% 'dirk34a', 'dirk22b', 'dirk33b')
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
    error('Invalid input argument : f must be an interger in {0,2,4,6,8,10} for Redonnet filter''s. (explicit)');
end

elseif strcmp(filtre,'shapiro') == 1
    ftr=1-sin(teta/2).^f;
    
elseif strcmp('visbal',filtre) == 1
    if f==0
        ftr=1;
    elseif f == 2
        num= (0.5+delta)+(0.5+delta)*cos(teta);
        denom=2*delta*cos(teta)+1;
        ftr=num./denom;
    elseif f == 4
        num= (5/8+3*delta/4)+(0.5+delta)*cos(teta)+(-1/8+delta/4)*cos(2*teta);
        denom=2*delta*cos(teta)+1;
        ftr=num./denom;
    elseif f == 6
        num= (11/16+5*delta/8)+(15/32+17*delta/16)*cos(teta)+(-3/16+3*delta/8)*cos(2*teta)+(1/32-delta/16)*cos(3*teta);
        denom=2*delta*cos(teta)+1;
        ftr=num./denom;
    elseif f == 8
        num= (93/128+70*delta/128)+(7/16+18*delta/16)*cos(teta)+(-7/32+14*delta/32)*cos(2*teta)+(1/16-delta/8)*cos(3*teta)+(-1/128+delta/64)*cos(4*teta);
        denom=2*delta*cos(teta)+1;
        ftr=num./denom;
    elseif f == 10
        num= (193/256+126*delta/256)+(105/256+302*delta/256)*cos(teta)+(15/64*(-1+2*delta))*cos(2*teta)+(45/512*(1-2*delta))*cos(3*teta)+(5/256*(-1+2*delta))*cos(4*teta)+(1-2*delta)/512*cos(5*teta);
        denom=2*delta*cos(teta)+1;
        ftr=num./denom;
    else
    error('Invalid input argument : f must be an interger in {0,2,4,6,8,10} for Visbal filter''s (implicit).');
    end
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
    s = 3*sin(teta)./(2+(cos(teta)));
elseif c == 6
    s = (4/3*sin(teta)-1/3*sin(2*teta))./(4/5+2*(2/15*cos(teta)-1/30*cos(2*teta)));
else
    error('Invalid input argument : c must be an integer in {4,6} for implicit space scheme.');
end

end

%% *** time discretisation (only with rk4) ********************************
% be careful : here i is the complex number such that 1i*1i = -1 !
if strcmp(time,'rk4') == 1
    g = 1-1i*lambda.*s-1/2*lambda.^2.*s.^2+1/6*1i*lambda.^3.*s.^3+1/24*lambda.^4.*s.^4;
    g=g.*ftr;
elseif strcmp(time,'rk1') == 1
    g=1-1i*lambda.*s;
    g=g.*ftr;
elseif strcmp(time,'rk2') == 1
    g = 1-1i*lambda.*s-1/2*lambda.^2.*s.^2;
    g=g.*ftr;
elseif strcmp(time,'dirk12a') == 1
    theta = -1i*lambda*s;
    g = (1+theta/2)./(1-theta/2);
    g = g.*ftr;
elseif strcmp(time,'dirk23a') == 1
    theta = -1i.*lambda.*s;
    nom = -2.*(18-9*theta-9*sqrt(3)*theta+3*theta.^3+2*sqrt(3).*theta.^3);
    denom = (-6 + 3.*theta + theta.*sqrt(3)).*(6-6*theta-2*sqrt(3).*theta+2*theta.^2+sqrt(3).*theta.^2);
    g=nom./denom;
    g = g.*ftr;
elseif strcmp(time,'dirk34a') == 1
    theta=-1i.*lambda.*s;
    alfa=2*cos(pi/18)/sqrt(3);
    num = -24*alfa+12*theta*alfa+36*theta*alfa.^2 +6.*theta.^2*alfa - 18*alfa^3*theta.^2 - 7*theta.^3.*alfa-9*theta.^3*alfa^2-3*theta.^3*alfa^3+3*theta.^3*alfa^4-2*theta.^3;
    denom = 3.*alfa.*(-2+theta+theta.*alfa).^3;
    g=num./denom;
    g=g.*ftr;
elseif strcmp(time,'dirk22b') == 1
    theta=-1i.*lambda.*s;
    ll=1-sqrt(2)/2;
    c2=1;
    num = 2*theta.^2.*ll^2-4*theta.*ll-4*theta.^2.*ll+2+theta.^2+2.*theta;
    denom = (-1+theta.*ll).^2;
    g=0.5*num./denom;
    g=g.*ftr;
elseif strcmp(time,'dirk33b') == 1
    theta=-1i.*lambda.*s;
    ll=0.4358665215;
    num=-8+24.*theta.*ll-42.*theta.^2.*ll.^2+6*theta.^2.*ll.^3-8.*theta+33.*theta.^2.*ll-5.*theta.^2;
    denom=(-1+theta.*ll).^3;
    g=1/8*num./denom;
    g=g.*ftr;
else
    error('Invalid input argument : time scheme (in time) is invalid');
end
    
    
end

