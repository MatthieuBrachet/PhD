function [fun,int] = fun_quad2(x,y,z,test)
global radius
global ang1 ang2
aaa=ang1; bbb=ang2;
x1=cos(aaa).*x-sin(aaa).*y;
y1=sin(aaa).*x+cos(aaa).*y;
z1=z;

x2=x1;
y2=cos(bbb).*y1-sin(bbb).*z1;
z2=sin(bbb).*y1+cos(bbb).*z1;

x=x2./radius; y=y2./radius; z=z2./radius;


if test == 0
    fun=ones(size(x));
    int=4*pi;
    
elseif test == 1
    fun=1+x+y.^2+y.*x.^2+x.^4+y.^5+(x.^2).*(y.^2).*(z.^2);
    int=216*pi/35;
    
elseif test == 2
    fun= 0.75*exp(-(9*x-2).^2/4-(9*y-2).^2/4-(9*z-2).^2/4)+...
        0.75*exp(-(9*x+1).^2/49-(9*y+1)/10-(9*z+1)/10)+...
        0.5*exp(-(9*x-7).^2/4-(9*y-3).^2/4-(9*z-5).^2/4) ...
        -0.2*exp(-(9*x-4).^2-(9*y-7).^2-(9*z-5).^2);
    int=6.6961822200736179523;
    
elseif test == 3
    fun=(1+tanh(-9*x-9*y+9*z))/9;
    int=4*pi/9;
    
elseif test == 4
    fun=(1+sign(-9*x-9*y+9*z))/9;
    int=4*pi/9;
    
elseif test == 5
    alpha=7;
    fun=(1-sign(pi*x+y))./alpha;
    int=4*pi/alpha;
end
int=int.*radius.^2;
end
