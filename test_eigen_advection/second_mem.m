function [ sm ] = second_mem( vect,time )
global n nn


[ funfI, funfII, funfIII, funfIV, funfV, funfVI ] = vect2fun( vect );

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

[vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time);

kI   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
kII  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
kIII = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
kIV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
kV   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
kVI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 

[ sm ] = fun2vect( kI, kII, kIII, kIV, kV, kVI );
end

