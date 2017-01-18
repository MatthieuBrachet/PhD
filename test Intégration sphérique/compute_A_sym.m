%%% COMPUTE THE MATRIX A SUCH THAT A*w=b WITH
%%% w BEING THE WEIGHTS FOR 1/8 PATCH (SYMMETRIES)
%%% err_i IS THE INTEGRALS OF SPHs FOR THE
%%% BASIC FORMULA (WITH THE METRIC TENSOR)

function [A,err_i] = compute_A_sym(nhs_max)
% nhs_max : nombre d'harmonique sphérique a intégrer correctement?

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
    if nhs_max==0,
        nhs_max=ceil(sqrt(6*N^2+2)-1);
        while mod(nhs_max,2)~=0,
            nhs_max=nhs_max+1;
        end      
    end

    %make_cs_grid(N);

    err_i=[];
    %%% construction of the matrix A
    A=[];
    for nhs=0:2:nhs_max,
        for mhs=0:4:nhs,
            a=zeros(1,nb_w);
            funfI=sph(nhs,mhs,x_fI,y_fI,z_fI);
            funfII=sph(nhs,mhs,x_fII,y_fII,z_fII);
            funfIII=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
            funfIV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
            funfV=sph(nhs,mhs,x_fV,y_fV,z_fV);
            funfVI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI);
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
        end
    end

end

