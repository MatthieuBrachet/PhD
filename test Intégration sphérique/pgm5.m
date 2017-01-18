% SHOW THE ERRORS OF THE NUMERICAL INTEGRATION ACCORDING TO nhs_max
% FOR f1, f2 AND f3

clear all;
global n nn;
global radius;
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=8;
nhs_max=2*N;
make_cs_grid(N);
    
[weights,opt_val]=eps_weights(nhs_max);

funfI=fun_f1(x_fI,y_fI,z_fI);
funfII=fun_f1(x_fII,y_fII,z_fII);
funfIII=fun_f1(x_fIII,y_fIII,z_fIII);
funfIV=fun_f1(x_fIV,y_fIV,z_fIV);
funfV=fun_f1(x_fV,y_fV,z_fV);
funfVI=fun_f1(x_fVI,y_fVI,z_fVI);
[nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
    int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);
abs(nrmg-216*pi/35)

funfI=fun_f2(x_fI,y_fI,z_fI);
funfII=fun_f2(x_fII,y_fII,z_fII);
funfIII=fun_f2(x_fIII,y_fIII,z_fIII);
funfIV=fun_f2(x_fIV,y_fIV,z_fIV);
funfV=fun_f2(x_fV,y_fV,z_fV);
funfVI=fun_f2(x_fVI,y_fVI,z_fVI);
[nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
    int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);
abs(nrmg-6.6961822200736179523)

funfI=fun_f3(x_fI,y_fI,z_fI);
funfII=fun_f3(x_fII,y_fII,z_fII);
funfIII=fun_f3(x_fIII,y_fIII,z_fIII);
funfIV=fun_f3(x_fIV,y_fIV,z_fIV);
funfV=fun_f3(x_fV,y_fV,z_fV);
funfVI=fun_f3(x_fVI,y_fVI,z_fVI);
[nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
    int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);
abs(nrmg-4*pi/9)    
    
    