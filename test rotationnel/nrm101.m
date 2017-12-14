function [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str)
global dxi deta dga;
global weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW VERSION COMPARED TO NRM72.M : COEFFT 1/4 FOR THE 8 CORNER POINTS
switch str
    case 'infty'
    nrmI=max(max(abs(funfI)));
    nrmII=max(max(abs(funfII)));
    nrmIII=max(max(abs(funfIII)));
    nrmIV=max(max(abs(funfIV)));
    nrmV=max(max(abs(funfV)));
    nrmVI=max(max(abs(funfVI)));
%
    nrmg=max([nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]);
    case '2'
        % face I
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*(funfI(2:nn-1,2:nn-1).^2)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*(funfI(1,2:nn-1).^2));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*(funfI(nn,2:nn-1).^2));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*(funfI(2:nn-1,1).^2));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*(funfI(2:nn-1,nn).^2));
    nrmNW=dxi*deta*dga(1,nn)*(funfI(1,nn)^2);
    nrmNE=dxi*deta*dga(nn,nn)*(funfI(nn,nn)^2);
    nrmSW=dxi*deta*dga(1,1)*(funfI(1,1)^2);
    nrmSE=dxi*deta*dga(nn,1)*(funfI(nn,1)^2);
    nrmI=sqrt(nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE));
        % face II
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*(funfII(2:nn-1,2:nn-1).^2)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*(funfII(1,2:nn-1).^2));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*(funfII(nn,2:nn-1).^2));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*(funfII(2:nn-1,1).^2));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*(funfII(2:nn-1,nn).^2));
    nrmNW=dxi*deta*dga(1,nn)*(funfII(1,nn)^2);
    nrmNE=dxi*deta*dga(nn,nn)*(funfII(nn,nn)^2);
    nrmSW=dxi*deta*dga(1,1)*(funfII(1,1)^2);
    nrmSE=dxi*deta*dga(nn,1)*(funfII(nn,1)^2);
    nrmII=sqrt(nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE));
            % face III
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*(funfIII(2:nn-1,2:nn-1).^2)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*(funfIII(1,2:nn-1).^2));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*(funfIII(nn,2:nn-1).^2));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*(funfIII(2:nn-1,1).^2));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*(funfIII(2:nn-1,nn).^2));
    nrmNW=dxi*deta*dga(1,nn)*(funfIII(1,nn)^2);
    nrmNE=dxi*deta*dga(nn,nn)*(funfIII(nn,nn)^2);
    nrmSW=dxi*deta*dga(1,1)*(funfIII(1,1)^2);
    nrmSE=dxi*deta*dga(nn,1)*(funfIII(nn,1)^2);
    nrmIII=sqrt(nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE));
                % face IV
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*(funfIV(2:nn-1,2:nn-1).^2)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*(funfIV(1,2:nn-1).^2));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*(funfIV(nn,2:nn-1).^2));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*(funfIV(2:nn-1,1).^2));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*(funfIV(2:nn-1,nn).^2));
    nrmNW=dxi*deta*dga(1,nn)*(funfIV(1,nn)^2);
    nrmNE=dxi*deta*dga(nn,nn)*(funfIV(nn,nn)^2);
    nrmSW=dxi*deta*dga(1,1)*(funfIV(1,1)^2);
    nrmSE=dxi*deta*dga(nn,1)*(funfIV(nn,1)^2);
    nrmIV=sqrt(nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE));
                   % face V
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*(funfV(2:nn-1,2:nn-1).^2)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*(funfV(1,2:nn-1).^2));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*(funfV(nn,2:nn-1).^2));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*(funfV(2:nn-1,1).^2));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*(funfV(2:nn-1,nn).^2));
    nrmNW=dxi*deta*dga(1,nn)*(funfV(1,nn)^2);
    nrmNE=dxi*deta*dga(nn,nn)*(funfV(nn,nn)^2);
    nrmSW=dxi*deta*dga(1,1)*(funfV(1,1)^2);
    nrmSE=dxi*deta*dga(nn,1)*(funfV(nn,1)^2);
    nrmV=sqrt(nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE)); 
                      % face VI
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*(funfV(2:nn-1,2:nn-1).^2)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*(funfVI(1,2:nn-1).^2));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*(funfVI(nn,2:nn-1).^2));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*(funfVI(2:nn-1,1).^2));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*(funfVI(2:nn-1,nn).^2));
    nrmNW=dxi*deta*dga(1,nn)*(funfVI(1,nn)^2);
    nrmNE=dxi*deta*dga(nn,nn)*(funfVI(nn,nn)^2);
    nrmSW=dxi*deta*dga(1,1)*(funfVI(1,1)^2);
    nrmSE=dxi*deta*dga(nn,1)*(funfVI(nn,1)^2);
    nrmVI=sqrt(nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE));  
    %
    nrmg=sqrt(nrmI^2+nrmII^2+nrmIII^2+nrmIV^2+nrmV^2+nrmVI^2);    
        case '1'
        % face I
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*abs(funfI(2:nn-1,2:nn-1))));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*abs(funfI(1,2:nn-1)));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*abs(funfI(nn,2:nn-1)));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*abs(funfI(2:nn-1,1)));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*abs(funfI(2:nn-1,nn)));
    nrmNW=dxi*deta*dga(1,nn)*abs(funfI(1,nn));
    nrmNE=dxi*deta*dga(nn,nn)*abs(funfI(nn,nn));
    nrmSW=dxi*deta*dga(1,1)*abs(funfI(1,1));
    nrmSE=dxi*deta*dga(nn,1)*abs(funfI(nn,1));
    nrmI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face II
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*abs(funfII(2:nn-1,2:nn-1))));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*abs(funfII(1,2:nn-1)));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*abs(funfII(nn,2:nn-1)));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*abs(funfII(2:nn-1,1)));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*abs(funfII(2:nn-1,nn)));
    nrmNW=dxi*deta*dga(1,nn)*abs(funfII(1,nn));
    nrmNE=dxi*deta*dga(nn,nn)*abs(funfII(nn,nn));
    nrmSW=dxi*deta*dga(1,1)*abs(funfII(1,1));
    nrmSE=dxi*deta*dga(nn,1)*abs(funfII(nn,1));
    nrmII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face III
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*abs(funfIII(2:nn-1,2:nn-1))));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*abs(funfIII(1,2:nn-1)));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*abs(funfIII(nn,2:nn-1)));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*abs(funfIII(2:nn-1,1)));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*abs(funfIII(2:nn-1,nn)));
    nrmNW=dxi*deta*dga(1,nn)*abs(funfIII(1,nn));
    nrmNE=dxi*deta*dga(nn,nn)*abs(funfIII(nn,nn));
    nrmSW=dxi*deta*dga(1,1)*abs(funfIII(1,1));
    nrmSE=dxi*deta*dga(nn,1)*abs(funfIII(nn,1));
    nrmIII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face IV
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*abs(funfIV(2:nn-1,2:nn-1))));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*abs(funfIV(1,2:nn-1)));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*abs(funfIV(nn,2:nn-1)));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*abs(funfIV(2:nn-1,1)));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*abs(funfIV(2:nn-1,nn)));
    nrmNW=dxi*deta*dga(1,nn)*abs(funfIV(1,nn));
    nrmNE=dxi*deta*dga(nn,nn)*abs(funfIV(nn,nn));
    nrmSW=dxi*deta*dga(1,1)*abs(funfIV(1,1));
    nrmSE=dxi*deta*dga(nn,1)*abs(funfIV(nn,1));
    nrmIV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face V
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*abs(funfV(2:nn-1,2:nn-1))));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*abs(funfV(1,2:nn-1)));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*abs(funfV(nn,2:nn-1)));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*abs(funfV(2:nn-1,1)));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*abs(funfV(2:nn-1,nn)));
    nrmNW=dxi*deta*dga(1,nn)*abs(funfV(1,nn));
    nrmNE=dxi*deta*dga(nn,nn)*abs(funfV(nn,nn));
    nrmSW=dxi*deta*dga(1,1)*abs(funfV(1,1));
    nrmSE=dxi*deta*dga(nn,1)*abs(funfV(nn,1));
    nrmV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face VI
    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*abs(funfVI(2:nn-1,2:nn-1))));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*abs(funfVI(1,2:nn-1)));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*abs(funfVI(nn,2:nn-1)));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*abs(funfVI(2:nn-1,1)));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*abs(funfVI(2:nn-1,nn)));
    nrmNW=dxi*deta*dga(1,nn)*abs(funfVI(1,nn));
    nrmNE=dxi*deta*dga(nn,nn)*abs(funfVI(nn,nn));
    nrmSW=dxi*deta*dga(1,1)*abs(funfVI(1,1));
    nrmSE=dxi*deta*dga(nn,1)*abs(funfVI(nn,1));
    nrmVI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
    %
    nrmg=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI;    
    
    case 'int' % INTEGRAL OF THE FUNCTION (IT IS NOT A NORM).
        % face I

    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*funfI(2:nn-1,2:nn-1)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*funfI(1,2:nn-1));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*funfI(nn,2:nn-1));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*funfI(2:nn-1,1));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*funfI(2:nn-1,nn));
    nrmNW=dxi*deta*dga(1,nn)*funfI(1,nn);
    nrmNE=dxi*deta*dga(nn,nn)*funfI(nn,nn);
    nrmSW=dxi*deta*dga(1,1)*funfI(1,1);
    nrmSE=dxi*deta*dga(nn,1)*funfI(nn,1);
    nrmI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face II

    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*funfII(2:nn-1,2:nn-1)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*funfII(1,2:nn-1));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*funfII(nn,2:nn-1));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*funfII(2:nn-1,1));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*funfII(2:nn-1,nn));
    nrmNW=dxi*deta*dga(1,nn)*funfII(1,nn);
    nrmNE=dxi*deta*dga(nn,nn)*funfII(nn,nn);
    nrmSW=dxi*deta*dga(1,1)*funfII(1,1);
    nrmSE=dxi*deta*dga(nn,1)*funfII(nn,1);
    nrmII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face III

    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*funfIII(2:nn-1,2:nn-1)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*funfIII(1,2:nn-1));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*funfIII(nn,2:nn-1));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*funfIII(2:nn-1,1));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*funfIII(2:nn-1,nn));
    nrmNW=dxi*deta*dga(1,nn)*funfIII(1,nn);
    nrmNE=dxi*deta*dga(nn,nn)*funfIII(nn,nn);
    nrmSW=dxi*deta*dga(1,1)*funfIII(1,1);
    nrmSE=dxi*deta*dga(nn,1)*funfIII(nn,1);
    nrmIII=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face IV

    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*funfIV(2:nn-1,2:nn-1)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*funfIV(1,2:nn-1));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*funfIV(nn,2:nn-1));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*funfIV(2:nn-1,1));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*funfIV(2:nn-1,nn));
    nrmNW=dxi*deta*dga(1,nn)*funfIV(1,nn);
    nrmNE=dxi*deta*dga(nn,nn)*funfIV(nn,nn);
    nrmSW=dxi*deta*dga(1,1)*funfIV(1,1);
    nrmSE=dxi*deta*dga(nn,1)*funfIV(nn,1);
    nrmIV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face V

    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*funfV(2:nn-1,2:nn-1)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*funfV(1,2:nn-1));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*funfV(nn,2:nn-1));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*funfV(2:nn-1,1));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*funfV(2:nn-1,nn));
    nrmNW=dxi*deta*dga(1,nn)*funfV(1,nn);
    nrmNE=dxi*deta*dga(nn,nn)*funfV(nn,nn);
    nrmSW=dxi*deta*dga(1,1)*funfV(1,1);
    nrmSE=dxi*deta*dga(nn,1)*funfV(nn,1);
    nrmV=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
            % face VI

    nrmint=dxi*deta*sum(sum(dga(2:nn-1,2:nn-1).*funfVI(2:nn-1,2:nn-1)));
    nrmW=dxi*deta*sum(dga(1,2:nn-1).*funfVI(1,2:nn-1));
    nrmE=dxi*deta*sum(dga(nn,2:nn-1).*funfVI(nn,2:nn-1));
    nrmS=dxi*deta*sum(dga(2:nn-1,1).*funfVI(2:nn-1,1));
    nrmN=dxi*deta*sum(dga(2:nn-1,nn).*funfVI(2:nn-1,nn));
    nrmNW=dxi*deta*dga(1,nn)*funfVI(1,nn);
    nrmNE=dxi*deta*dga(nn,nn)*funfVI(nn,nn);
    nrmSW=dxi*deta*dga(1,1)*funfVI(1,1);
    nrmSE=dxi*deta*dga(nn,1)*funfVI(nn,1);
    nrmVI=nrmint+(1/2)*(nrmW+nrmE+nrmS+nrmN)+(1/3)*(nrmNW+nrmNE+nrmSW+nrmSE);
    %
    nrmg=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI; 
    
    case 'cor_int' % INTEGRALE CORRIGEE DE BRICE PORTENELLE
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
    
    case 'simpson' % INTRAGE CONSERVATIVE BASEE SUR SIMPSON
        
    if mod(n+1,2)==1
        error('n+1 must be odd for this quadrature formula');
    end
     
    A=8*ones(nn,nn); 
    A(2:2:end-1,2:2:end-1)=16;
    A(1:2:end,1:2:end)=4;

    A(1,1:2:end)=2;
    A(end,1:2:end)=2;
    A(1:2:end,1)=2;
    A(1:2:end,end)=2;

    A(1,2:2:end-1)=4;
    A(end,2:2:end-1)=4;
    A(2:2:end-1,1)=4;
    A(2:2:end-1,end)=4;

    A(1,1)=1;
    A(1,end)=1;
    A(end,1)=1;
    A(end,end)=1;

    wei=A./9;
    
    nrmI=1^2.*dxi.*deta.*sum(sum(dga.*wei.*funfI));
    nrmII=1^2.*dxi.*deta.*sum(sum(dga.*wei.*funfII));
    nrmIII=1^2.*dxi.*deta.*sum(sum(dga.*wei.*funfIII));
    nrmIV=1^2.*dxi.*deta.*sum(sum(dga.*wei.*funfIV));
    nrmV=1^2.*dxi.*deta.*sum(sum(dga.*wei.*funfV));
    nrmVI=1^2.*dxi.*deta.*sum(sum(dga.*wei.*funfVI));
    nrmg=nrmI+nrmII+nrmIII+nrmIV+nrmV+nrmVI; 
        
end









