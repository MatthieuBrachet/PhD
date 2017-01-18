function [A,err_i] = compute_A(nhs_max)
%%% COMPUTE THE MATRIX A SUCH THAT A*w=b WITH
%%% w BEING THE WEIGHTS FOR 1 PATCH
%%% err_i IS THE INTEGRALS OF SPHs FOR THE
%%% BASIC FORMULA (WITH THE METRIC TENSOR)

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
    if nhs_max==0,
        nhs_max=ceil(sqrt(6*N^2+2)-1);
        while mod(nhs_max,2)~=0,
            nhs_max=nhs_max+1;
        end      
    end

    err_i=[];
    A=[];
    for nhs=0:nhs_max,
        for mhs=0:nhs, 
            a=zeros(1,nn*nn);
            funfI=sph(nhs,mhs,x_fI,y_fI,z_fI);
            funfII=sph(nhs,mhs,x_fII,y_fII,z_fII);
            funfIII=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
            funfIV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
            funfV=sph(nhs,mhs,x_fV,y_fV,z_fV);
            funfVI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI);
            err_i=[err_i;int_weights(dxi*deta*dga,funfI,funfII,funfIII,funfIV,funfV,funfVI)];
            funfI=reshape(funfI',[nn*nn,1]);
            funfII=reshape(funfII',[nn*nn,1]);
            funfIII=reshape(funfIII',[nn*nn,1]);
            funfIV=reshape(funfIV',[nn*nn,1]);
            funfV=reshape(funfV',[nn*nn,1]);
            funfVI=reshape(funfVI',[nn*nn,1]);
            for j=1:nn*nn,
                if j==1 || j==nn || j==(nn-1)*nn+1 || j==nn*nn,
                    coeff=1/3;
                elseif (j>1 && j<nn) || (j>(nn-1)*nn+1 && j<nn*nn) || ...
                            (j>1 && j<(nn-1)*nn+1 && mod(j-1,nn)==0) || ...
                            (j>nn && j<nn*nn && mod(j,nn)==0),
                    coeff=1/2;
                else
                    coeff=1;
                end 
                a(j)=coeff*(funfI(j)+funfII(j)+funfIII(j)+funfIV(j)+funfV(j)+funfVI(j));
            end
            A=[A;a];
        end
    end

end

