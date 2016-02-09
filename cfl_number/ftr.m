function [ ftr ] = ftr( teta,f,filtre,delta )
%% ************************************************************************
% author : Brachet Matthieu (IECL-Metz)
% birth date : janv-26-2016
% last version : 08-feb-2016
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

