function [ scaf1f2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI ] = ...
    intscaC_weights( weights,fun1fI,fun1fII,fun1fIII,fun1fIV,fun1fV,fun1fVI,...
    fun2fI,fun2fII,fun2fIII,fun2fIV,fun2fV,fun2fVI)
% ROUTINE DEDUCED FROM INT_WEIGHTS
    nn=size(weights,1);
    % SALAR PRODUCT OF THE FUNCTION F1 AND F2. 
        % face I
    %disp('fI');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*fun1fI(2:nn-1,2:nn-1).*conj(fun2fI(2:nn-1,2:nn-1))));
    nrmW=sum(weights(1,2:nn-1).*fun1fI(1,2:nn-1).*conj(fun2fI(1,2:nn-1)));
    nrmE=sum(weights(nn,2:nn-1).*fun1fI(nn,2:nn-1).*conj(fun2fI(nn,2:nn-1)));
    nrmS=sum(weights(2:nn-1,1).*fun1fI(2:nn-1,1).*conj(fun2fI(2:nn-1,1)));
    nrmN=sum(weights(2:nn-1,nn).*fun1fI(2:nn-1,nn).*conj(fun2fI(2:nn-1,nn)));
    nrmNW=weights(1,nn)*fun1fI(1,nn)*conj(fun2fI(1,nn));
    nrmNE=weights(nn,nn)*fun1fI(nn,nn)*conj(fun2fI(nn,nn));
    nrmSW=weights(1,1)*fun1fI(1,1)*conj(fun2fI(1,1));
    nrmSE=weights(nn,1)*fun1fI(nn,1)*conj(fun2fI(nn,1));
    nrmI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face II
    %disp('fII');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*fun1fII(2:nn-1,2:nn-1).*conj(fun2fII(2:nn-1,2:nn-1))));
    nrmW=sum(weights(1,2:nn-1).*fun1fII(1,2:nn-1).*conj(fun2fII(1,2:nn-1)));
    nrmE=sum(weights(nn,2:nn-1).*fun1fII(nn,2:nn-1).*conj(fun2fII(nn,2:nn-1)));
    nrmS=sum(weights(2:nn-1,1).*fun1fII(2:nn-1,1).*conj(fun2fII(2:nn-1,1)));
    nrmN=sum(weights(2:nn-1,nn).*fun1fII(2:nn-1,nn).*conj(fun2fII(2:nn-1,nn)));
    nrmNW=weights(1,nn)*fun1fII(1,nn)*conj(fun2fII(1,nn));
    nrmNE=weights(nn,nn)*fun1fII(nn,nn)*conj(fun2fII(nn,nn));
    nrmSW=weights(1,1)*fun1fII(1,1)*conj(fun2fII(1,1));
    nrmSE=weights(nn,1)*fun1fII(nn,1)*conj(fun2fII(nn,1));
    nrmII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);  
         % face III
    %disp('fIII');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*fun1fIII(2:nn-1,2:nn-1).*conj(fun2fIII(2:nn-1,2:nn-1))));
    nrmW=sum(weights(1,2:nn-1).*fun1fIII(1,2:nn-1).*conj(fun2fIII(1,2:nn-1)));
    nrmE=sum(weights(nn,2:nn-1).*fun1fIII(nn,2:nn-1).*conj(fun2fIII(nn,2:nn-1)));
    nrmS=sum(weights(2:nn-1,1).*fun1fIII(2:nn-1,1).*conj(fun2fIII(2:nn-1,1)));
    nrmN=sum(weights(2:nn-1,nn).*fun1fIII(2:nn-1,nn).*conj(fun2fIII(2:nn-1,nn)));
    nrmNW=weights(1,nn)*fun1fIII(1,nn)*conj(fun2fIII(1,nn));
    nrmNE=weights(nn,nn)*fun1fIII(nn,nn)*conj(fun2fIII(nn,nn));
    nrmSW=weights(1,1)*fun1fIII(1,1)*conj(fun2fIII(1,1));
    nrmSE=weights(nn,1)*fun1fIII(nn,1)*conj(fun2fIII(nn,1));
    nrmIII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);   
            % face IV
    %disp('fIV');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*fun1fIV(2:nn-1,2:nn-1).*conj(fun2fIV(2:nn-1,2:nn-1))));
    nrmW=sum(weights(1,2:nn-1).*fun1fIV(1,2:nn-1).*conj(fun2fIV(1,2:nn-1)));
    nrmE=sum(weights(nn,2:nn-1).*fun1fIV(nn,2:nn-1).*conj(fun2fIV(nn,2:nn-1)));
    nrmS=sum(weights(2:nn-1,1).*fun1fIV(2:nn-1,1).*conj(fun2fIV(2:nn-1,1)));
    nrmN=sum(weights(2:nn-1,nn).*fun1fIV(2:nn-1,nn).*conj(fun2fIV(2:nn-1,nn)));
    nrmNW=weights(1,nn)*fun1fIV(1,nn)*conj(fun2fIV(1,nn));
    nrmNE=weights(nn,nn)*fun1fIV(nn,nn)*conj(fun2fIV(nn,nn));
    nrmSW=weights(1,1)*fun1fIV(1,1)*conj(fun2fIV(1,1));
    nrmSE=weights(nn,1)*fun1fIV(nn,1)*conj(fun2fIV(nn,1));
    nrmIV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
           % face V
    %disp('fV');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*fun1fV(2:nn-1,2:nn-1).*conj(fun2fV(2:nn-1,2:nn-1))));
    nrmW=sum(weights(1,2:nn-1).*fun1fV(1,2:nn-1).*conj(fun2fV(1,2:nn-1)));
    nrmE=sum(weights(nn,2:nn-1).*fun1fV(nn,2:nn-1).*conj(fun2fV(nn,2:nn-1)));
    nrmS=sum(weights(2:nn-1,1).*fun1fV(2:nn-1,1).*conj(fun2fV(2:nn-1,1)));
    nrmN=sum(weights(2:nn-1,nn).*fun1fV(2:nn-1,nn).*conj(fun2fV(2:nn-1,nn)));
    nrmNW=weights(1,nn)*fun1fV(1,nn)*conj(fun2fV(1,nn));
    nrmNE=weights(nn,nn)*fun1fV(nn,nn)*conj(fun2fV(nn,nn));
    nrmSW=weights(1,1)*fun1fV(1,1)*conj(fun2fV(1,1));
    nrmSE=weights(nn,1)*fun1fV(nn,1)*conj(fun2fV(nn,1));
    nrmV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
           % face VI
    %disp('fVI');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*fun1fVI(2:nn-1,2:nn-1).*conj(fun2fVI(2:nn-1,2:nn-1))));
    nrmW=sum(weights(1,2:nn-1).*fun1fVI(1,2:nn-1).*conj(fun2fVI(1,2:nn-1)));
    nrmE=sum(weights(nn,2:nn-1).*fun1fVI(nn,2:nn-1).*conj(fun2fVI(nn,2:nn-1)));
    nrmS=sum(weights(2:nn-1,1).*fun1fVI(2:nn-1,1).*conj(fun2fVI(2:nn-1,1)));
    nrmN=sum(weights(2:nn-1,nn).*fun1fVI(2:nn-1,nn).*conj(fun2fVI(2:nn-1,nn)));
    nrmNW=weights(1,nn)*fun1fVI(1,nn)*conj(fun2fVI(1,nn));
    nrmNE=weights(nn,nn)*fun1fVI(nn,nn)*conj(fun2fVI(nn,nn));
    nrmSW=weights(1,1)*fun1fVI(1,1)*conj(fun2fVI(1,1));
    nrmSE=weights(nn,1)*fun1fVI(nn,1)*conj(fun2fVI(nn,1));
    nrmVI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
    %
    % RESULT IS THE SCALAR PRODUCT OF FUNCTIONS F1 AND F2
    scaf1f2=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI;    
end

