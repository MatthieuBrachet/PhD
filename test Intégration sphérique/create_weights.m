function weights = create_weights( n,choice )
    
    weights = zeros(n,n);
    
    if choice==1,
        % CREATION OF THE WEIGHTS
        % with these properties : 
        % 1) m(i,-j)=m(i,j)
        % 3) m(j,i)=m(i,j)
        % these 2 properties => 2) m(-i,j)=m(i,j)
        % m(i,j)>0
        for i=1:(n-1)/2+1,
            for j=1:i,
                weights(i,j)=abs(randn(1));
                weights(j,i)=weights(i,j);
                weights(i,n-j+1)=weights(i,j);
                weights(j,n-i+1)=weights(j,i);
                
                weights(n-j+1,i)=weights(i,n-j+1);
                weights(n-i+1,j)=weights(j,n-i+1);
                weights(n-j+1,n-i+1)=weights(n-j+1,i);
                weights(n-i+1,n-j+1)=weights(n-i+1,j);
            end
        end  
        
        
    elseif choice==2,
        % CREATION OF THE WEIGHTS
        % with these properties : 
        % 1) m(i,-j)=m(i,j)
        % 2) m(-i,j)=m(i,j)
        % m(i,j)>0
        for i=1:(n-1)/2+1,
            for j=1:(n-1)/2+1,
                weights(i,j)=abs(randn(1));
                weights(i,n-j+1)=weights(i,j);
                weights(n-i+1,j)=weights(i,j);
                weights(n-i+1,n-j+1)=weights(i,j);
            end
        end  
        
    elseif choice==3,
        % CREATION OF THE WEIGHTS
        % with these properties : 
        % 3) m(j,i)=m(i,j)
        % m(i,j)>0
        % => SYMMETRIC MATRICE
        for i=1:n,
            for j=1:i,
                weights(i,j)=abs(randn(1));
                weights(j,i)=weights(i,j);
            end
        end  
        
    elseif choice==4,
        % CREATION OF THE WEIGHTS
        % with these properties : 
        % m(n+1-i+1,n+1-j+1)=m(i,j)
        % m(i,j)>0
        % => PERSYMMETRIC MATRICE
        for i=1:n,
            for j=n-i+1:-1:1,
                weights(i,j)=abs(randn(1));
            end
        end 
        weights=rot90(weights);
        for i=1:n,
            for j=1:i-1,
                weights(j,i)=weights(i,j);
            end
        end 
        weights=rot90(weights,-1);
        
    elseif choice==5,
        % CREATION OF THE WEIGHTS
        % with these properties : 
        % m(j,i)=m(i,j)
        % m(n+1-i+1,n+1-j+1)=m(i,j)
        % m(i,j)>0
        % => SP MATRICE (symmetric and persymmetric)
        for i=1:n/2+1,
            for j=1:i,
                weights(i,j)=abs(randn(1));
                weights(j,i)=weights(i,j);
                weights(n+1-i+1,n+1-j+1)=weights(i,j);
                weights(n+1-j+1,n+1-i+1)=weights(i,j);
            end
        end
        for i=n/2+2:n+1,
            for j=1:n/2+1-(i-(n/2+1)),
                weights(i,j)=abs(randn(1));
                weights(j,i)=weights(i,j);
                weights(n+1-i+1,n+1-j+1)=weights(i,j);
                weights(n+1-j+1,n+1-i+1)=weights(i,j);
            end
        end
        
        
end

