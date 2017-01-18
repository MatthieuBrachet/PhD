%%% CALCULUS OF THE INTEGRAL OVER S2 FOR DIFFERENT FUNCTIONS %%%
%%% WARNING : IT IS RELATIVE ERROR !!!

function err = int_funs_fornberg( weights )

    global n nn;
    global x_fI y_fI z_fI;
    global x_fII y_fII z_fII;
    global x_fIII y_fIII z_fIII;
    global x_fIV y_fIV z_fIV;
    global x_fV y_fV z_fV;
    global x_fVI y_fVI z_fVI;
    global dxi deta dga;
    global radius
    
    w=weights;
    err=[];

    % % FUN_CST: =1 everywhere
    %
    funfI=fun_cst(x_fI,y_fI,z_fI);
    funfII=fun_cst(x_fII,y_fII,z_fII);
    funfIII=fun_cst(x_fIII,y_fIII,z_fIII);
    funfIV=fun_cst(x_fIV,y_fIV,z_fIV);
    funfV=fun_cst(x_fV,y_fV,z_fV);
    funfVI=fun_cst(x_fVI,y_fVI,z_fVI);
    [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
       int_weights(w,funfI,funfII,funfIII,funfIV,funfV,funfVI);
    %disp 'erreur fonction constante = 1'
    err=[err,abs(nrmg-4*pi*radius.^2)./abs(4*pi*radius.^2)];

    % % FUN_F1: f1 in Fornberg-Martel p.1176
    %
    funfI=fun_f1(x_fI,y_fI,z_fI);
    funfII=fun_f1(x_fII,y_fII,z_fII);
    funfIII=fun_f1(x_fIII,y_fIII,z_fIII);
    funfIV=fun_f1(x_fIV,y_fIV,z_fIV);
    funfV=fun_f1(x_fV,y_fV,z_fV);
    funfVI=fun_f1(x_fVI,y_fVI,z_fVI);
    [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
       int_weights(w,funfI,funfII,funfIII,funfIV,funfV,funfVI);
    %disp 'erreur f1'
    err=[err,abs(nrmg-radius.^2*216*pi/35)./abs(radius.^2*216*pi/35)];

%     % % FUN_F2: f2 in Fornberg-Martel p.1176
%     %
%     funfI=fun_f2(x_fI,y_fI,z_fI);
%     funfII=fun_f2(x_fII,y_fII,z_fII);
%     funfIII=fun_f2(x_fIII,y_fIII,z_fIII);
%     funfIV=fun_f2(x_fIV,y_fIV,z_fIV);
%     funfV=fun_f2(x_fV,y_fV,z_fV);
%     funfVI=fun_f2(x_fVI,y_fVI,z_fVI);
%     [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
%        int_weights(w,funfI,funfII,funfIII,funfIV,funfV,funfVI);
%     %disp 'erreur f2'
%     err=[err,abs(nrmg-6.6961822200736179523*radius.^2)./abs(6.6961822200736179523*radius.^2)];

    % % FUN_F3: f3 in Fornberg-Martel p.1176
    %
    funfI=fun_f3(x_fI,y_fI,z_fI);
    funfII=fun_f3(x_fII,y_fII,z_fII);
    funfIII=fun_f3(x_fIII,y_fIII,z_fIII);
    funfIV=fun_f3(x_fIV,y_fIV,z_fIV);
    funfV=fun_f3(x_fV,y_fV,z_fV);
    funfVI=fun_f3(x_fVI,y_fVI,z_fVI);
    [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
       int_weights(w,funfI,funfII,funfIII,funfIV,funfV,funfVI);
    %disp 'erreur f3'
    err=[err,abs(nrmg-4*pi/9*radius.^2)./abs(4*pi/9*radius.^2)];

    % % FUN_F4: f4 in Fornberg-Martel p.1176
    %
    funfI=fun_f4(x_fI,y_fI,z_fI);
    funfII=fun_f4(x_fII,y_fII,z_fII);
    funfIII=fun_f4(x_fIII,y_fIII,z_fIII);
    funfIV=fun_f4(x_fIV,y_fIV,z_fIV);
    funfV=fun_f4(x_fV,y_fV,z_fV);
    funfVI=fun_f4(x_fVI,y_fVI,z_fVI);
    [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
       int_weights(w,funfI,funfII,funfIII,funfIV,funfV,funfVI);
    %disp 'erreur f4'
    err=[err,abs(nrmg-4*pi/9*radius.^2)./abs(4*pi/9*radius.^2)];

end

