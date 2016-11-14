function [ V_xi, V_eta ] = riemann74( h, u_xi, u_eta )
global gp

G11=dot(u_xi,u_xi);
G12=dot(u_xi,u_eta);
G22=dot(u_eta,u_eta);

V11=-G12;
V12=gp;
V13=gp;
V21=G11;
V22=0;
V23=0;
V31=0;
V32=sqrt(G11*gp*h);
V33=-sqrt(G11*gp*h);
V_xi=[V11 V12 V13; V21 V22 V23; V31 V32 V33];

V11=G22;
V12=0;
V13=0;
V21=-G12;
V22=sqrt(G22*h*gp);
V23=-sqrt(G22*h*gp);
V31=0;
V32=gp;
V33=gp;
V_eta=[V11 V12 V13; V21 V22 V23; V31 V32 V33];


end

