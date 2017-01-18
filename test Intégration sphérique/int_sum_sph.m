% SUM OF THE NUMERICAL INTEGRALS (abs) OF SPH
% FOR nhs EVEN UNTIL nhs_max AND mhs<=nhs,mhs=0[4]

function res = int_sum_sph( weights,nhs_max )

    global x_fI y_fI z_fI;
    global x_fII y_fII z_fII;
    global x_fIII y_fIII z_fIII;
    global x_fIV y_fIV z_fIV;
    global x_fV y_fV z_fV;
    global x_fVI y_fVI z_fVI;

    res = 0;
    
    for nhs=2:2:nhs_max,
        for mhs=0:4:nhs,            
            funfI=sph(nhs,mhs,x_fI,y_fI,z_fI);
            funfII=sph(nhs,mhs,x_fII,y_fII,z_fII);
            funfIII=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
            funfIV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
            funfV=sph(nhs,mhs,x_fV,y_fV,z_fV);
            funfVI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI); 
            res = res + abs(int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI));
        end
    end    

end

