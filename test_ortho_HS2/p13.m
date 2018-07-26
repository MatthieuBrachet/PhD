clc; clear all; close all;
format long;
N=2;
n=N-1;nn=N+1;
make_cs_grid(N);
eps=0.1;

harm=[];
dim=[];
for nhs=0:20
    for mhs=-nhs:nhs
        [ hs ] = sph_cs( nhs, mhs );
        harm=[harm hs];
    end
    dim=[dim rank(harm,eps)];
end
[1 diff(dim)]


nhs=nhs+1;
disp(['nhs = ' num2str(nhs)]) 
[ hs ] = sph_cs( nhs, 0 );
harm=[harm hs];
dim=[dim rank(harm,eps)];
disp('ajout de |mhs| = 0')
disp(['gain : ' num2str(dim(end)-dim(end-1))])

disp('***************************')
for mhs=1:nhs
    [ hs ] = sph_cs( nhs, mhs );
    harm=[harm hs];
    [ hs ] = sph_cs( nhs, -mhs );
    harm=[harm hs];
    dim=[dim rank(harm,eps)];
    disp(['ajout de |mhs| = ' num2str(mhs)])
    disp(['gain : ' num2str(dim(end)-dim(end-1))])
    disp('***************************')
end