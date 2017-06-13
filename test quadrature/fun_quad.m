function [fun,int] = fun_quad(x,y,z,test)
global radius dga
x=x./radius; y=y./radius; z=z./radius;

if test == -1
    nhs=10;
    mhs=8;
    fun=sph(nhs,mhs,x,y,z)+eps;
    int=eps*4*pi;

elseif test == 0
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
    int=0.414319590372970; % for panel I
    % int=6.6961822200736179523;
    
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
    
elseif test == 6
    fun=radius^2./dga;
    %int=2.467401100272339; % PANEL I
    int=14.804406601634035;
elseif test == 7
    fun=(radius^2./dga).^4;
    %int=17.579738281473187; %^2
    int= 25.368886759078510;
end
int=int.*radius.^2;
end

