function [grad_fI,grad_fII,grad_fIII,grad_fIV,grad_fV,grad_fVI]=...
    grad1000(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
global pgxi kgxi
global pgeta kgeta
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;

%% --- PANEL I
Mfunf=reshape(mfunfI,[],1);
du=pgxi\(kgxi*Mfunf);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*Mfunf);
dueta=reshape(du,nn,nn);

grad_fI(:,:,1)=duxi.*gxi_I(:,:,1)+dueta.*geta_I(:,:,1);
grad_fI(:,:,2)=duxi.*gxi_I(:,:,2)+dueta.*geta_I(:,:,2);
grad_fI(:,:,3)=duxi.*gxi_I(:,:,3)+dueta.*geta_I(:,:,3);

%% --- PANEL II
Mfunf=reshape(mfunfII,[],1);
du=pgxi\(kgxi*Mfunf);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*Mfunf);
dueta=reshape(du,nn,nn);

grad_fII(:,:,1)=duxi.*gxi_II(:,:,1)+dueta.*geta_II(:,:,1);
grad_fII(:,:,2)=duxi.*gxi_II(:,:,2)+dueta.*geta_II(:,:,2);
grad_fII(:,:,3)=duxi.*gxi_II(:,:,3)+dueta.*geta_II(:,:,3);

%% --- PANEL III
Mfunf=reshape(mfunfIII,[],1);
du=pgxi\(kgxi*Mfunf);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*Mfunf);
dueta=reshape(du,nn,nn);

grad_fIII(:,:,1)=duxi.*gxi_III(:,:,1)+dueta.*geta_III(:,:,1);
grad_fIII(:,:,2)=duxi.*gxi_III(:,:,2)+dueta.*geta_III(:,:,2);
grad_fIII(:,:,3)=duxi.*gxi_III(:,:,3)+dueta.*geta_III(:,:,3);

%% --- PANEL IV
Mfunf=reshape(mfunfIV,[],1);
du=pgxi\(kgxi*Mfunf);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*Mfunf);
dueta=reshape(du,nn,nn);

grad_fIV(:,:,1)=duxi.*gxi_IV(:,:,1)+dueta.*geta_IV(:,:,1);
grad_fIV(:,:,2)=duxi.*gxi_IV(:,:,2)+dueta.*geta_IV(:,:,2);
grad_fIV(:,:,3)=duxi.*gxi_IV(:,:,3)+dueta.*geta_IV(:,:,3);

%% --- PANEL V
Mfunf=reshape(mfunfV,[],1);
du=pgxi\(kgxi*Mfunf);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*Mfunf);
dueta=reshape(du,nn,nn);

grad_fV(:,:,1)=duxi.*gxi_V(:,:,1)+dueta.*geta_V(:,:,1);
grad_fV(:,:,2)=duxi.*gxi_V(:,:,2)+dueta.*geta_V(:,:,2);
grad_fV(:,:,3)=duxi.*gxi_V(:,:,3)+dueta.*geta_V(:,:,3);
    
%% --- PANEL VI
Mfunf=reshape(mfunfVI,[],1);
du=pgxi\(kgxi*Mfunf);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*Mfunf);
dueta=reshape(du,nn,nn);

grad_fVI(:,:,1)=duxi.*gxi_VI(:,:,1)+dueta.*geta_VI(:,:,1);
grad_fVI(:,:,2)=duxi.*gxi_VI(:,:,2)+dueta.*geta_VI(:,:,2);
grad_fVI(:,:,3)=duxi.*gxi_VI(:,:,3)+dueta.*geta_VI(:,:,3);





end

