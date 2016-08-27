function[adv_fI,adv_fII,adv_fIII, adv_fIV, adv_fV, adv_fVI]=...
    adv72(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)


for i=1:3
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=gr72(funfI(:,:,i),funfII(:,:,i),funfIII(:,:,i),funfIV(:,:,i),funfV(:,:,i),funfVI(:,:,i),n,nn);
    
    adv_fI(:,:,i)=funfI(:,:,1).*grad_I(:,:,1)+funfI(:,:,2).*grad_I(:,:,2)+funfI(:,:,3).*grad_I(:,:,3);
    adv_fII(:,:,i)=funfII(:,:,1).*grad_II(:,:,1)+funfII(:,:,2).*grad_II(:,:,2)+funfII(:,:,3).*grad_II(:,:,3);
    adv_fIII(:,:,i)=funfIII(:,:,1).*grad_III(:,:,1)+funfI(:,:,2).*grad_III(:,:,2)+funfIII(:,:,3).*grad_III(:,:,3);
    adv_fIV(:,:,i)=funfIV(:,:,1).*grad_IV(:,:,1)+funfIV(:,:,2).*grad_IV(:,:,2)+funfIV(:,:,3).*grad_IV(:,:,3);
    adv_fV(:,:,i)=funfV(:,:,1).*grad_V(:,:,1)+funfV(:,:,2).*grad_V(:,:,2)+funfV(:,:,3).*grad_V(:,:,3);
    adv_fVI(:,:,i)=funfVI(:,:,1).*grad_VI(:,:,1)+funfVI(:,:,2).*grad_VI(:,:,2)+funfVI(:,:,3).*grad_VI(:,:,3);
    
end

end