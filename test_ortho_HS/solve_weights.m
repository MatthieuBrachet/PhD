%%% Solve the linear system A(1:k,:)*w=b, 
%%% with err=0 if we want the total weights, and
%%% err~=0 if we want the epsilon weights ;
%%% sym=0 if we want all the weights for one patch, and
%%% sym~=0 if we want only weights for 1/8 of the patch, and then
%%% reconstruct the total weights by symmetries

%%% If A==0, it constructs the matrix A
%%% If k==0, it chooses k=min(size(A))

function res = solve_weights( A,err_i,k,err,sym )

    global n nn;
    global radius;
    global dxi deta dga;
    global x_fI y_fI z_fI;
    global x_fII y_fII z_fII;
    global x_fIII y_fIII z_fIII;
    global x_fIV y_fIV z_fIV;
    global x_fV y_fV z_fV;
    global x_fVI y_fVI z_fVI;

    % Computation of the matrix A according to sym
    if A==0,
        if sym==0,
            [A,err_i]=compute_A();
        else
            [A,err_i]=compute_A_sym();
        end
    end    
    
    % Computation of the vector b according to err
    if k==0,
        %k=min(size(A));
        k=size(A,1);
    end
    b=zeros(k,1);
    if err==0,
        b(1)=sqrt(4*pi);
    else
        b(1)=4*pi*radius^2/sqrt(4*pi)-err_i(1);
        for i=2:k,
            b(i)=-err_i(i);
        end
    end

    % Computation of tmp, the solution of A(1:k,:)*tmp=b
    tmp=pinv(A(1:k,:))*b;
    
    % Reshape of the weights
    N=nn-1;
    if sym==0,
        w=reshape(tmp,[nn,nn])';
    else
        w=zeros(nn,nn);
        ii=1;
        for i=1:N/2,
            for j=1:i,
                % .><.  
                w(i,j)=tmp(ii);
                w(i,nn-j+1)=w(i,j);
                w(nn-i+1,j)=w(i,j);
                w(nn-i+1,nn-j+1)=w(i,j);        
                % >:<
                w(j,i)=w(i,j);
                w(j,nn-i+1)=w(j,i);
                w(nn-j+1,i)=w(j,i);
                w(nn-j+1,nn-i+1)=w(j,i);

                ii=ii+1;
            end
        end   
        for i=1:N/2,
            % middle line and middle column, 
            % except the middle element of the matrix
            w((nn-1)/2+1,i)=tmp(ii);
            w((nn-1)/2+1,nn-i+1)=w((nn-1)/2+1,i);
            w(i,(nn-1)/2+1)=w((nn-1)/2+1,i);
            w(nn-i+1,(nn-1)/2+1)=w(i,(nn-1)/2+1);
            ii=ii+1;
        end
        % middle element of the matrix
        w(N/2+1,N/2+1)=tmp(ii);
    end

    % Result, according to if we want w or epsilon_w
    if err==0,
        res=w;
    else
        res=w./dxi./deta;
    end

end

