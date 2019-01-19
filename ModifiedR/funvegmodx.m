function Ut=funvegmodx(t,U)


global kappa Q f Sb Hd Vw Nb 
global Gamma Lev wC  phiZr Db M Kb
global Nx Dx1 Dx2 
global P Tyear alpha Cnrm

tmp=reshape(U,Nx,3);

H=tmp(:,1);
S=tmp(:,2);
B=tmp(:,3);


Jx= -Vw*(H.^Hd)./(1.0+Nb*B);
In= kappa*(B + Q*f)./(B + Q).*(H).*(1.0-S).^Sb;

Ht= -(Dx1*Jx) -  In +  precip(P,t,Tyear,alpha,Cnrm); 
St=  (-(Lev+Gamma*B).*S + In)/phiZr;
Bt=    Db*(Dx2*B ) +   (wC*Gamma*S.*(1.0-B/Kb)-M).*B ;

%display(precip(p,t))
Ut=[Ht(:); St(:); Bt(:)];

function pt=precip(p,t,Tyear,alpha,Cnrm)
%pt=p  *  sech(  cos(pi*t/Tyear) / epsilon ).^2  * Constant
 pt=p  *  sech(alpha*cos(pi*t/Tyear)           ).^2  *  Cnrm;
 %pt=p  *  sech(8*cos(pi*t/Tyear)           ).^2  *  12.482485478454578897;
 
%these values may only be accurate to 3 or 4 decimals.
%epsilon=1, Constant=1.494515826022178
%epsilon=1/2, Constant=2.806475267441896
%epsilon=1/4, Constant=6.088198059309244
%epsilon=1/8, Constant=12.481052357761302
%epsilon=1/16, Constant=25.085971409699642

