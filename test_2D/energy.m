function [ E ] = energy(h,u,v)
global gp hp
E=sum(.5*gp.*h.^2+.5*hp.*(u.^2+v.^2));
end

