function [ pdt ] = scal(u,v)
pdt=u(:,:,1).*v(:,:,1)+u(:,:,2).*v(:,:,2)+u(:,:,3).*v(:,:,3);
end