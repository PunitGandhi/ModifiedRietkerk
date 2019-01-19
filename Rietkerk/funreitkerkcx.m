function Ut=funvegmodx(t,U)


global Dh Dw  
global a f p g m v
global Nx
global Dx1 Dx2 
global p0 Tyear

tmp=reshape(U,Nx,3);

H=tmp(:,1);
W=tmp(:,2);
B=tmp(:,3);

Gw= W./(1.0+W);
In= a*(B + f)./(B + 1.0);

Ht= Dh*(Dx2*H) -  In.*H + precip(t,p0,Tyear);  

Wt= Dw*(Dx2*W) + In.*H - v*W - g*Gw.*B;

Bt= Dx2*B - m*B + Gw.*B;

%display(precip(p,t))
Ut=[Ht(:); Wt(:); Bt(:)];

function pt=precip(t,p0,Tyear)
%pt=p0  *  sech(  cos(pi*t/Tyear) / epsilon ).^2  * Constant
 pt=p0  *  sech(8*cos(pi*t/Tyear)           ).^2  *  12.482485478454578897;
 
%these values may only be accurate to 3 or 4 decimals.
%epsilon=1, Constant=1.494515826022178
%epsilon=1/2, Constant=2.806475267441896
%epsilon=1/4, Constant=6.088198059309244
%epsilon=1/8, Constant=12.481052357761302
%epsilon=1/16, Constant=25.085971409699642

