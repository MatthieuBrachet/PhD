close all
%% courbes

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
        div72( vte_fI  , vte_fII  , vte_fIII  , vte_fIV  , vte_fV  , vte_fVI  , n , nn);
    
figure(8)
plot_cs11(n,nn,div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI);
%%
[divb_fI,divb_fII,divb_fIII,divb_fIV,divb_fV,divb_fVI]=...
        div72b( vte_fI  , vte_fII  , vte_fIII  , vte_fIV  , vte_fV  , vte_fVI  , n , nn);
    
figure(9)
plot_cs11(n,nn,divb_fI,divb_fII,divb_fIII,divb_fIV,divb_fV,divb_fVI);
%%
figure(10);
plot_cs11(n,nn,divb_fI-div_fI,divb_fII-div_fII,divb_fIII-div_fIII,divb_fIV-div_fIV,divb_fV-div_fV,divb_fVI-div_fVI);