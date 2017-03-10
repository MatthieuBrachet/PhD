clc; clear all; close all;
% fun1 is spherical harmonic arounf (Oz).
global dxi deta dga;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

N=4;
make_cs_grid(N);

% Computation of weights with basic formula
weights=dxi*deta*dga;

hs=6;
for nhs=1:length(hs)
    nhs_max=hs(nhs);
    NHS=[0:1:nhs_max];
    for i=1:length(NHS)
        clc; [NHS(i) nhs_max]
        nhs1=NHS(i);
        MHS1=[-nhs1:1:nhs1];
        for p=1:length(MHS1)
            mhs1=MHS1(p);

            fun1fI=conj(sph(nhs1,mhs1,x_fI,y_fI,z_fI));
            fun1fII=conj(sph(nhs1,mhs1,x_fII,y_fII,z_fII));
            fun1fIII=conj(sph(nhs1,mhs1,x_fIII,y_fIII,z_fIII));
            fun1fIV=conj(sph(nhs1,mhs1,x_fIV,y_fIV,z_fIV));
            fun1fV=conj(sph(nhs1,mhs1,x_fV,y_fV,z_fV));
            fun1fVI=conj(sph(nhs1,mhs1,x_fVI,y_fVI,z_fVI));


            for j=1:length(NHS)
                nhs2=NHS(j);
                MHS2=[-nhs2:1:nhs2];
                for q=1:length(MHS2)
                    mhs2=MHS2(q);

                    fun2fI=sph(nhs2,mhs2,x_fI,y_fI,z_fI);
                    fun2fII=sph(nhs2,mhs2,x_fII,y_fII,z_fII);
                    fun2fIII=sph(nhs2,mhs2,x_fIII,y_fIII,z_fIII);
                    fun2fIV=sph(nhs2,mhs2,x_fIV,y_fIV,z_fIV);
                    fun2fV=sph(nhs2,mhs2,x_fV,y_fV,z_fV);
                    fun2fVI=sph(nhs2,mhs2,x_fVI,y_fVI,z_fVI);

                    [scaf1f2,nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI]=...
                        intsca_weights( weights,fun1fI,fun1fII,fun1fIII,fun1fIV,fun1fV,fun1fVI,...
                        fun2fI,fun2fII,fun2fIII,fun2fIV,fun2fV,fun2fVI);

                    mat_scalar(i+p-1,j+q-1)=scaf1f2;

                end
            end
        end
    end

    rg(nhs)=rank(mat_scalar);
end

figure(1)
plot(hs,rg,hs,(6*N.^2+2).*ones(size(hs)))
legend('rank','max')
xlabel('nhs_max')
ylabel('rank')
