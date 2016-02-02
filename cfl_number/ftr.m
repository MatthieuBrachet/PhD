function [ ftr ] = ftr( teta,f,filtre,delta )
%% ************************************************************************
% author : Brachet Matthieu (IECL-Metz)
% birth date : janv-26-2016
% last version : janv-27-2016
%
%% ************************************************************************
% this function build the amplification number obtain with use of a linear
% filter.
%
%% *** inputs
% teta : wave number;
% f : order of the space filter;
% filtre : type of space filter.
%
%% *** outputs
% ftr : amplification number.
%
%% ************************************************************************

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

end

