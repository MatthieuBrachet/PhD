function [ nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI ] = ...
    int_weights( weights,funfI,funfII,funfIII,funfIV,funfV,funfVI )

    nn=size(weights,1);

    % INTEGRAL OF THE FUNCTION (IT IS NOT A NORM).
        % face I
    %disp('fI');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*funfI(2:nn-1,2:nn-1)));
    nrmW=sum(weights(1,2:nn-1).*funfI(1,2:nn-1));
    nrmE=sum(weights(nn,2:nn-1).*funfI(nn,2:nn-1));
    nrmS=sum(weights(2:nn-1,1).*funfI(2:nn-1,1));
    nrmN=sum(weights(2:nn-1,nn).*funfI(2:nn-1,nn));
    nrmNW=weights(1,nn)*funfI(1,nn);
    nrmNE=weights(nn,nn)*funfI(nn,nn);
    nrmSW=weights(1,1)*funfI(1,1);
    nrmSE=weights(nn,1)*funfI(nn,1);
    nrmI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face II
    %disp('fII');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*funfII(2:nn-1,2:nn-1)));
    nrmW=sum(weights(1,2:nn-1).*funfII(1,2:nn-1));
    nrmE=sum(weights(nn,2:nn-1).*funfII(nn,2:nn-1));
    nrmS=sum(weights(2:nn-1,1).*funfII(2:nn-1,1));
    nrmN=sum(weights(2:nn-1,nn).*funfII(2:nn-1,nn));
    nrmNW=weights(1,nn)*funfII(1,nn);
    nrmNE=weights(nn,nn)*funfII(nn,nn);
    nrmSW=weights(1,1)*funfII(1,1);
    nrmSE=weights(nn,1)*funfII(nn,1);
    nrmII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face III
    %disp('fIII');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*funfIII(2:nn-1,2:nn-1)));
    nrmW=sum(weights(1,2:nn-1).*funfIII(1,2:nn-1));
    nrmE=sum(weights(nn,2:nn-1).*funfIII(nn,2:nn-1));
    nrmS=sum(weights(2:nn-1,1).*funfIII(2:nn-1,1));
    nrmN=sum(weights(2:nn-1,nn).*funfIII(2:nn-1,nn));
    nrmNW=weights(1,nn)*funfIII(1,nn);
    nrmNE=weights(nn,nn)*funfIII(nn,nn);
    nrmSW=weights(1,1)*funfIII(1,1);
    nrmSE=weights(nn,1)*funfIII(nn,1);
    nrmIII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face IV
    %disp('fIV');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*funfIV(2:nn-1,2:nn-1)));
    nrmW=sum(weights(1,2:nn-1).*funfIV(1,2:nn-1));
    nrmE=sum(weights(nn,2:nn-1).*funfIV(nn,2:nn-1));
    nrmS=sum(weights(2:nn-1,1).*funfIV(2:nn-1,1));
    nrmN=sum(weights(2:nn-1,nn).*funfIV(2:nn-1,nn));
    nrmNW=weights(1,nn)*funfIV(1,nn);
    nrmNE=weights(nn,nn)*funfIV(nn,nn);
    nrmSW=weights(1,1)*funfIV(1,1);
    nrmSE=weights(nn,1)*funfIV(nn,1);
    nrmIV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face V
    %disp('fV');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*funfV(2:nn-1,2:nn-1)));
    nrmW=sum(weights(1,2:nn-1).*funfV(1,2:nn-1));
    nrmE=sum(weights(nn,2:nn-1).*funfV(nn,2:nn-1));
    nrmS=sum(weights(2:nn-1,1).*funfV(2:nn-1,1));
    nrmN=sum(weights(2:nn-1,nn).*funfV(2:nn-1,nn));
    nrmNW=weights(1,nn)*funfV(1,nn);
    nrmNE=weights(nn,nn)*funfV(nn,nn);
    nrmSW=weights(1,1)*funfV(1,1);
    nrmSE=weights(nn,1)*funfV(nn,1);
    nrmV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face VI
    %disp('fVI');
    nrmint=sum(sum(weights(2:nn-1,2:nn-1).*funfVI(2:nn-1,2:nn-1)));
    nrmW=sum(weights(1,2:nn-1).*funfVI(1,2:nn-1));
    nrmE=sum(weights(nn,2:nn-1).*funfVI(nn,2:nn-1));
    nrmS=sum(weights(2:nn-1,1).*funfVI(2:nn-1,1));
    nrmN=sum(weights(2:nn-1,nn).*funfVI(2:nn-1,nn));
    nrmNW=weights(1,nn)*funfVI(1,nn);
    nrmNE=weights(nn,nn)*funfVI(nn,nn);
    nrmSW=weights(1,1)*funfVI(1,1);
    nrmSE=weights(nn,1)*funfVI(nn,1);
    nrmVI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
    %
    nrmg=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI;    
end

