function [weights,cvx_optval] = eps_weights( nhs_max )

    global n nn;
    global radius;
    global dxi deta dga;
    global x_fI y_fI z_fI;
    global x_fII y_fII z_fII;
    global x_fIII y_fIII z_fIII;
    global x_fIV y_fIV z_fIV;
    global x_fV y_fV z_fV;
    global x_fVI y_fVI z_fVI;
    
    funfI=fun_cst(x_fI,y_fI,z_fI);
    funfII=fun_cst(x_fII,y_fII,z_fII);
    funfIII=fun_cst(x_fIII,y_fIII,z_fIII);
    funfIV=fun_cst(x_fIV,y_fIV,z_fIV);
    funfV=fun_cst(x_fV,y_fV,z_fV);
    funfVI=fun_cst(x_fVI,y_fVI,z_fVI);

    cvx_begin
        variable weights(nn,nn) symmetric
        minimize int_sum_sph(weights,nhs_max)
        subject to
            for i=1:(nn-1)/2,
                % 4 interiors of the matrix
                for j=1:i,
                    % .><.
                    weights(i,nn-j+1)==weights(i,j)
                    weights(nn-i+1,j)==weights(i,j)
                    weights(nn-i+1,nn-j+1)==weights(i,j)     
                    % >:<
                    weights(j,nn-i+1)==weights(j,i)
                    weights(nn-j+1,i)==weights(j,i)
                    weights(nn-j+1,nn-i+1)==weights(j,i)
                end
                % middle line and middle column, 
                % except the middle element of the matrix
                weights((nn-1)/2+1,nn-i+1)==weights((nn-1)/2+1,i)
                weights(i,(nn-1)/2+1)==weights((nn-1)/2+1,i)
                weights(nn-i+1,(nn-1)/2+1)==weights(i,(nn-1)/2+1)
            end
            weights >= 0
            int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI) == 4*pi*radius^2                  
    cvx_end


end

