function [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_mixte74(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
kappa=1;

%% filtrage d'ordre 2n=8
% ftr = 8
[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_beta74a(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
[funftI8,funftII8,funftIII8,funftIV8,funftV8,funftVI8]=...
    ftr_alpha74a(funftI,funftII,funftIII,funftIV,funftV,funftVI,n,nn);

%% perturbation produit mixte
% ftr1 = 4
[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_beta74b(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
u_fI   = funftI  - funfI;
u_fII  = funftII - funfII;
u_fIII = funftIII- funfIII;
u_fIV  = funftIV - funfIV;
u_fV   = funftV  - funfV;
u_fVI  = funftVI - funfVI;

[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_alpha74b(u_fI,u_fII,u_fIII,u_fIV,u_fV,u_fVI,n,nn);
uf_fI   = funftI   - u_fI;
uf_fII  = funftII  - u_fII;
uf_fIII = funftIII - u_fIII;
uf_fIV  = funftIV  - u_fIV;
uf_fV   = funftV   - u_fV;
uf_fVI  = funftVI  - u_fVI;

%% assemblage
funftI   = funftI8   - kappa*uf_fI;
funftII  = funftII8  - kappa*uf_fII;
funftIII = funftIII8 - kappa*uf_fIII;
funftIV  = funftIV8  - kappa*uf_fIV;
funftV   = funftV8   - kappa*uf_fV;
funftVI  = funftVI8  - kappa*uf_fVI;

% funftI=funftI8;
% funftII=funftII8;
% funftIII=funftIII8;
% funftIV=funftIV8;
% funftV=funftV8;
% funftVI=funftVI8;

end

