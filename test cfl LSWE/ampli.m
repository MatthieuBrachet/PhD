function [ stab ] = ampli(teta1, teta2,dx,dy,dt)



%% constantes
% coriolis
teta0=pi/4;
omega = 7.27*10^-5;
f0 = 2*omega*cos(teta0);

% athmosphere
H=1;

% gravity
g=9.81;


for i=1:length(teta1)
    for j=1:length(teta2)
        t1=teta1(i);
        t2=teta2(j);
        
        %% discretisations
        s1 = sin(t1)./(2/3+2*(1/6*cos(t1)));
        s2 = sin(t2)./(2/3+2*(1/6*cos(t2)));

        %% filtre
        ftr1 = 772/1024+2*(210/1024*cos(t1)-120/1024*cos(2*t1)+45/1024*cos(3*t1)-10/1024*cos(4*t1)+1/1024*cos(5*t1));
        ftr2 = 772/1024+2*(210/1024*cos(t2)-120/1024*cos(2*t2)+45/1024*cos(3*t2)-10/1024*cos(4*t2)+1/1024*cos(5*t2));

        %% matrice
        A=[0  f0  -1i*g/dx*s1; -f0  0  -1i*g/dy*s2; -1i*H/dx*s1  -1i*H/dy*s2  0];

        %% coefficient d'amplification
        G=eye(3,3)+dt*A+1/2*dt^2*A^2+1/6*dt^3*A^3+1/24*dt^4*A^4;


        %% calcul du parametre de stabilité
        spectre=eigs(G)*ftr1*ftr2;
        stab(i,j)=max(abs(spectre));
    end
end

