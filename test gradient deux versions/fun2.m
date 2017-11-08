function [func,gr] = fun2(x,y,z)
global radius
p=1;
q=2;
r=3;
func=(x.^p).*(y.^q).*(z.^r);
proj=(p+q+r)./radius.*func;

gr(:,:,1)=p.*x.^(p-1).*y.^q.*z.^r-proj.*x./radius;
gr(:,:,2)=q.*x.^p.*y.^(q-1).*z.^r-proj.*y./radius;
gr(:,:,3)=r.*x.^p.*y.^q.*z.^(r-1)-proj.*z./radius;
end

