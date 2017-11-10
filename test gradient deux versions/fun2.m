function [func,gr] = fun2(x,y,z)
global radius
rr=radius;
p=1;
q=1;
r=1;
func=(x.^p).*(y.^q).*(z.^r);
proj=(p+q+r)./rr.*func;

gr(:,:,1)=p.*x.^(p-1).*y.^q.*z.^r-proj.*x./rr;
gr(:,:,2)=q.*x.^p.*y.^(q-1).*z.^r-proj.*y./rr;
gr(:,:,3)=r.*x.^p.*y.^q.*z.^(r-1)-proj.*z./rr;
end

