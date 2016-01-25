function [ g ] = ark4cf(lambda, teta,c,space,f,filtre,alpha)
%% ************************************************************************
% author : Brachet Matthieu (IECL-Metz)
% birth date : janv-19-2016
% last version : janv-19-2016
%
%% ************************************************************************
% this function build the amplification number obtain with rk4 in space
% discretisation and an other scheme for space discretisation and a space
% filter on the advection equation.
%
%% *** name of this function
% rk4 : time discretisation rk4;
% c : space discretisation;
% f0 :pas de filtre.
%
%% *** inputs
% lambda : cfl number;
% teta : wave number;
% c : order for the space scheme.
% space : type of discretisation (explicit or implicit).
% f : order of the space filter;
% filtre : type of space filter (filter = 0 : explicit, filter = 1,
% implicit).
%
%% *** outputs
% g : amplification number.
%
%% ************************************************************************

%% *** calcul du filtre ***************************************************
if strcmp(filtre,'explicit') == 1
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

elseif strcmp(filtre,'implicit') == 1
if f==0
    ftr=1;
elseif f==2
    ftr = 0.5+2*(1/4*cos(teta));
    ftr=ftr./(1+2*alpha*cos(teta));
elseif f==4
    ftr = 10/16+2*(4/16*cos(teta)-1/16*cos(2*teta));
    ftr=ftr./(1+2*alpha*cos(teta));
elseif f==6
    ftr = 44/64+2*(15/64*cos(teta)-6/64*cos(2*teta)+1/64*cos(3*teta));
    ftr=ftr./(1+2*alpha*cos(teta));
elseif f == 8
    ftr = 186/256+2*(56/256*cos(teta)-28/256*cos(2*teta)+8/256*cos(3*teta)-1/256*cos(4*teta));
    ftr=ftr./(1+2*alpha*cos(teta));
elseif f == 10
    ftr = 772/1024+2*(210/1024*cos(teta)-120/1024*cos(2*teta)+45/1024*cos(3*teta)-10/1024*cos(4*teta)+1/1024*cos(5*teta));
    ftr=ftr./(1+2*alpha*cos(teta));
else
    error('Invalid input argument : f must be an interger in {0,2,4,6,8,10} for explicit filter.');
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
    s = sin(teta)./(2/3+2*(1/6*cos(teta)));
elseif c == 6
    s = (4/3*sin(teta)-1/3*sin(2*teta))./(4/5+2*(2/15*cos(teta)-1/30*cos(2*teta)));
else
    error('Invalid input argument : c must be an integer in {4,6} for implicit space scheme.');
end

end

%% *** time discretisation (only with rk4) ********************************
% be careful : here i is the complex number such that i*i = -1 !
g = 1-i*lambda.*s-1/2*lambda.^2.*s.^2+1/6*i*lambda.^3.*s.^3+1/24*lambda.^4.*s.^4;
g=g.*ftr;
end

