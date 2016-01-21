function [ g ] = ark4cf(lambda, teta,c,space,f,filtre)
%% ************************************************************************
% author :
%       - Brachet Matthieu (IECL-Metz)
% birth date : janv-19-2016
% last version : janv-19-2016
%
%% ************************************************************************
% this function build the amplification number obtain with rk4 in space
% discretisation and an other scheme for space discretisation and a space
% filter on the advection equation.
%% *** name of this function
% a : amplification coefficient;
% rk4 : time discretisation rk4;
% c : space discretisation;
% f0 :pas de filtre.
%
%% *** inputs
% lambda : cfl number;
% teta : wave number;
% c : order for the space scheme.
% space : type of discretisation (space = 0 : explicit, space = 1,
% implicit)
% f : order of the space filter;
% filtre : type of space filter (space = 0 : explicit, space = 1,
% implicit)
%
%% *** outputs
% g : amplification number.


%% *** calcul du filtre ***************************************************
if f==0
    ftr=1;
elseif f==2



end

