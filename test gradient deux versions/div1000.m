function [div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    grad1000(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn)
global pgxi kgxi
global pgeta kgeta

%% --- PANEL I
MfunfI=reshape(mfunfI,[],1);
du=pgxi\(kgxi*MfunfI);
duxi=reshape(du,nn,nn);

MfunfI=reshape(mfunfI,[],1);
du=pgxi\(kgxi*MfunfI);
duxi=reshape(du,nn,nn);
du=pgeta\(kgeta*MfunfI);
dueta=reshape(du,nn,nn);


    
    






end

