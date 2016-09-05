function [adv_fI,adv_fII,adv_fIII,adv_fIV,adv_fV,adv_fVI]=...
    adv72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
%% NE FONCTIONNE PAS !!!

%% *** terme 1 ************************************************************
%% *** RESEAU 1 ***
[n1,n2,~]=size(mfunfI);
for i=1:n1
    for j=1:n2
        u=[mfunfI(i,j,1) mfunfI(i,j,2) mfunfI(i,j,3)];
        tens_I(i,j,:,:)=u'*u;
    end
end

%% *** RESEAU 2 ***
[n1,n2,~]=size(mfunfII);
for i=1:n1
    for j=1:n2
        u=[mfunfII(i,j,1) mfunfII(i,j,2) mfunfII(i,j,3)];
        tens_II(i,j,:,:)=u'*u;
    end
end

%% *** RESEAU 3 ***
[n1,n2,~]=size(mfunfIII);
for i=1:n1
    for j=1:n2
        u=[mfunfIII(i,j,1) mfunfIII(i,j,2) mfunfIII(i,j,3)];
        tens_III(i,j,:,:)=u'*u;
    end
end

%% *** RESEAU 4 ***
[n1,n2,~]=size(mfunfIV);
for i=1:n1
    for j=1:n2
        u=[mfunfIV(i,j,1) mfunfIV(i,j,2) mfunfIV(i,j,3)];
        tens_IV(i,j,:,:)=u'*u;
    end
end

%% *** RESEAU 5 ***
[n1,n2,~]=size(mfunfV);
for i=1:n1
    for j=1:n2
        u=[mfunfV(i,j,1) mfunfV(i,j,2) mfunfV(i,j,3)];
        tens_V(i,j,:,:)=u'*u;
    end
end

%% *** RESEAU 6 ***
[n1,n2,~]=size(mfunfVI);
for i=1:n1
    for j=1:n2
        u=[mfunfVI(i,j,1) mfunfVI(i,j,2) mfunfVI(i,j,3)];
        tens_VI(i,j,:,:)=u'*u;
    end
end

%% *** divergence du tenseur

for k=1:3
    [ter_1I(:,:,k),ter_1II(:,:,k),ter_1III(:,:,k),ter_1IV(:,:,k),ter_1V(:,:,k),ter_1VI(:,:,k)]=div72(tens_I(:,:,k,:),tens_II(:,:,k,:), tens_III(:,:,k,:), tens_IV(:,:,k,:), tens_V(:,:,k,:), tens_VI(:,:,k,:),n,nn);
end


%% *** terme 2 ************************************************************
[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div72(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

for i=1:n1
    for j=1:n2
        for k=1:3
            ter_2I(i,j,k)=div_fI(i,j)*mfunfI(i,j,k);
            ter_2II(i,j,k)=div_fII(i,j)*mfunfII(i,j,k);
            ter_2III(i,j,k)=div_fIII(i,j)*mfunfIII(i,j,k);
            ter_2IV(i,j,k)=div_fIV(i,j)*mfunfIV(i,j,k);
            ter_2V(i,j,k)=div_fV(i,j)*mfunfV(i,j,k);
            ter_2VI(i,j,k)=div_fVI(i,j)*mfunfVI(i,j,k);
        end
    end
end

%% *** ASSEMBLAGE *********************************************************

adv_fI=ter_1I-ter_2I;
adv_fII=ter_1II-ter_2II;
adv_fIII=ter_1III-ter_2III;
adv_fIV=ter_1IV-ter_2IV;
adv_fV=ter_1V-ter_2V;
adv_fVI=ter_1VI-ter_2VI;


end

