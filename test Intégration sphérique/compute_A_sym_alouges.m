%%% COMPUTE THE MATRIX A SUCH THAT A*w=b WITH
%%% w BEING THE WEIGHTS FOR 1/8 PATCH (SYMMETRIES)
%%% err_i IS THE INTEGRALS OF SPHs FOR THE
%%% BASIC FORMULA (WITH THE METRIC TENSOR)

function [A,err_i] = compute_A_sym_alouges(r,nt,nl)

    global n nn;
    global radius;
    global dxi deta dga;
    global x_fI y_fI z_fI;
    global x_fII y_fII z_fII;
    global x_fIII y_fIII z_fIII;
    global x_fIV y_fIV z_fIV;
    global x_fV y_fV z_fV;
    global x_fVI y_fVI z_fVI;

    N=nn-1;
    nb_w=N^2/8+3*N/4+1;
   
    make_cs_grid(N);

    err_i=[];
    %%% construction of the matrix A
    A=[];
    
    %%%%%%%%%%%%%%%%%%
%     funfI=fun_cst(x_fI,y_fI,z_fI);
%     funfII=fun_cst(x_fII,y_fII,z_fII);
%     funfIII=fun_cst(x_fIII,y_fIII,z_fIII);
%     funfIV=fun_cst(x_fIV,y_fIV,z_fIV);
%     funfV=fun_cst(x_fV,y_fV,z_fV);
%     funfVI=fun_cst(x_fVI,y_fVI,z_fVI);
%     err_i=[err_i;int_weights(dxi*deta*dga,funfI,funfII,funfIII,funfIV,funfV,funfVI)];
%     tmp=funfI+funfII+funfIII+funfIV+funfV+funfVI;
%     j=1;
%     for ii=1:N/2+1,
%         for jj=1:ii,
%             if ii==1,
%                 coeff=1/3;
%             elseif (ii>1 && jj==1),
%                 coeff=1/2;
%             else
%                 coeff=1;
%             end
%             % diagonales (sauf point central)
%             if (ii==jj && ii~=N/2+1),
%                 a(j)=coeff*(tmp(ii,jj)+tmp(ii,nn-jj+1)+tmp(nn-ii+1,jj)+tmp(nn-ii+1,nn-jj+1));
%             % abscisses et ordonnées (sauf point central)
%             elseif (ii==N/2+1 && ii~=jj),
%                 a(j)=coeff*(tmp(ii,jj)+tmp(jj,ii)+tmp(ii,nn-jj+1)+tmp(nn-jj+1,ii));
%             % autres cas du triangle (1/8 de panel) (sauf point
%             % central)
%             elseif (ii~=N/2+1 || jj~=N/2+1),
%                 a(j)=coeff*(tmp(ii,jj)+tmp(nn-ii+1,jj)+tmp(ii,nn-jj+1)+tmp(nn-ii+1,nn-jj+1)+...
%                         tmp(jj,ii)+tmp(jj,nn-ii+1)+tmp(nn-jj+1,ii)+tmp(nn-jj+1,nn-ii+1));
%             % point central
%             else 
%                 a(j)=coeff*(tmp(ii,jj));
%             end
%             j=j+1;                
%         end
%     end
%     A=[A;a];
%     %%%%%%%%%%%%%%%%%%%
    
    theta=0;
    lambda=0;
    %rr=r;
    for rr=0:r,
    %for theta=linspace(-pi/2,pi/2,nt),
        %for lambda=linspace(pi,pi,nl),
            x=rr*cos(theta)*cos(lambda);
            y=rr*cos(theta)*sin(lambda);
            z=rr*sin(theta);
            X=[x;y;z];
            a=zeros(1,nb_w);
            funfI=fun_alouges(X,x_fI,y_fI,z_fI);
            funfII=fun_alouges(X,x_fII,y_fII,z_fII);
            funfIII=fun_alouges(X,x_fIII,y_fIII,z_fIII);
            funfIV=fun_alouges(X,x_fIV,y_fIV,z_fIV);
            funfV=fun_alouges(X,x_fV,y_fV,z_fV);
            funfVI=fun_alouges(X,x_fVI,y_fVI,z_fVI);
            err_i=[err_i;int_weights(dxi*deta*dga,funfI,funfII,funfIII,funfIV,funfV,funfVI)];
            tmp=funfI+funfII+funfIII+funfIV+funfV+funfVI;
            j=1;
            for ii=1:N/2+1,
                for jj=1:ii,
                    if ii==1,
                        coeff=1/3;
                    elseif (ii>1 && jj==1),
                        coeff=1/2;
                    else
                        coeff=1;
                    end
                    % diagonales (sauf point central)
                    if (ii==jj && ii~=N/2+1),
                        a(j)=coeff*(tmp(ii,jj)+tmp(ii,nn-jj+1)+tmp(nn-ii+1,jj)+tmp(nn-ii+1,nn-jj+1));
                    % abscisses et ordonnées (sauf point central)
                    elseif (ii==N/2+1 && ii~=jj),
                        a(j)=coeff*(tmp(ii,jj)+tmp(jj,ii)+tmp(ii,nn-jj+1)+tmp(nn-jj+1,ii));
                    % autres cas du triangle (1/8 de panel) (sauf point
                    % central)
                    elseif (ii~=N/2+1 || jj~=N/2+1),
                        a(j)=coeff*(tmp(ii,jj)+tmp(nn-ii+1,jj)+tmp(ii,nn-jj+1)+tmp(nn-ii+1,nn-jj+1)+...
                                tmp(jj,ii)+tmp(jj,nn-ii+1)+tmp(nn-jj+1,ii)+tmp(nn-jj+1,nn-ii+1));
                    % point central
                    else 
                        a(j)=coeff*(tmp(ii,jj));
                    end
                    j=j+1;                
                end
            end
            A=[A;a];
        %end
    %end
    end

end