%%% View of the integrals of all the SPH of order (nhs,mhs) with 0<=mhs<=nhs
%%% on the sphere ;
%%% If detail~=0, then it shows the integrals for the 6 patches
%%% If sph4~=0, then it only shows the SPH of order (nhs,mhs) with
%%% mhs = 0 mod 4

function err_i = view_int_sph( nhs,weights,detail,sph4 )

    teta0=0;
    lambda0=0;

    global x_fI y_fI z_fI;
    global x_fII y_fII z_fII;
    global x_fIII y_fIII z_fIII;
    global x_fIV y_fIV z_fIV;
    global x_fV y_fV z_fV;
    global x_fVI y_fVI z_fVI;

    err_i=[];
    err_fI=[];err_fII=[];err_fIII=[];
    err_fIV=[];err_fV=[];err_fVI=[];
    ii_a=[];
    step=1;
    if sph4~=0,
        step=4;
    end
    for mhs=-nhs:step:nhs,
%         funfI=sph(nhs,mhs,x_fI,y_fI,z_fI);
%         funfII=sph(nhs,mhs,x_fII,y_fII,z_fII);
%         funfIII=sph(nhs,mhs,x_fIII,y_fIII,z_fIII);
%         funfIV=sph(nhs,mhs,x_fIV,y_fIV,z_fIV);
%         funfV=sph(nhs,mhs,x_fV,y_fV,z_fV);
%         funfVI=sph(nhs,mhs,x_fVI,y_fVI,z_fVI);
        funfI=sph_rot(nhs,mhs,x_fI,y_fI,z_fI,teta0,lambda0);
        funfII=sph_rot(nhs,mhs,x_fII,y_fII,z_fII,teta0,lambda0);
        funfIII=sph_rot(nhs,mhs,x_fIII,y_fIII,z_fIII,teta0,lambda0);
        funfIV=sph_rot(nhs,mhs,x_fIV,y_fIV,z_fIV,teta0,lambda0);
        funfV=sph_rot(nhs,mhs,x_fV,y_fV,z_fV,teta0,lambda0);
        funfVI=sph_rot(nhs,mhs,x_fVI,y_fVI,z_fVI,teta0,lambda0);

        % APPROXIMATE QUADRATURE ON THE SPHERE
        [nrmg,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
            int_weights(weights,funfI,funfII,funfIII,funfIV,funfV,funfVI);
        ii_a=[ii_a,nrmg];
        ii_ex=0; % EXACT INTEGRAL OF SPHERICAL HARMONICS
        err_i=[err_i,abs(nrmg-ii_ex)];
        err_fI=[err_fI,nrmI];
        err_fII=[err_fII,nrmII];
        err_fIII=[err_fIII,nrmIII];
        err_fIV=[err_fIV,nrmIV];
        err_fV=[err_fV,nrmV];
        err_fVI=[err_fVI,nrmVI];
    end
    if detail == 1
        err_i
    elseif detail == 2
        err_i
        err_fI
        err_fII
        err_fIII
        err_fIV
        err_fV
        err_fVI
    end


end

