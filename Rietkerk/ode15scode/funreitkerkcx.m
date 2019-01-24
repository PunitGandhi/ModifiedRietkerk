function Ut=funreitkerkcx(t,U,parameters)

Dh = parameters.Dh; % diffusion of surface water
Dw = parameters.Dw; % diffusion of soil water.
a = parameters.a;  % infiltration parameter
f = parameters.f;  % infiltration parameter
g = parameters.g;  % transpiration param
v = parameters.v;  % evaporation (nu)
m = parameters.m;  % mortality
Dx1 = parameters.Dx1;  % differentiation matrix for first derivative - not used here
Dx2 = parameters.Dx2;  % diff. matrix for second derivative
p0 = parameters.p0;    % mean annual precipitation
Tyear = parameters.Tyear;  % length of a year in units of equations
Nx = parameters.Nx;  % number of spatial discretization points
pseason = parameters.pseason;  % precip parameter - seasonality
Cnrm = parameters.Cnrm;        % normalization constant


tmp=reshape(U,Nx,3);

H=tmp(:,1);
W=tmp(:,2);
B=tmp(:,3);

Gw= W./(1.0+W);
In= a*(B + f)./(B + 1.0);

Ht= Dh*(Dx2*H) -  In.*H + precip(t,p0,Tyear, pseason, Cnrm);  

Wt= Dw*(Dx2*W) + In.*H - v*W - g*Gw.*B;

Bt= Dx2*B - m*B + Gw.*B;

%display(precip(p,t))
Ut=[Ht(:); Wt(:); Bt(:)];

function pt=precip(t,p0,Tyear, pseason, Cnrm)
%pt=p0  *  sech(  cos(pi*t/Tyear) / epsilon ).^2  * Constant
 pt=p0  *  sech(pseason*cos(pi*t/Tyear)           ).^2  *  Cnrm;
 
%these values may only be accurate to 3 or 4 decimals.
%epsilon=1, Constant=1.494515826022178
%epsilon=1/2, Constant=2.806475267441896
%epsilon=1/4, Constant=6.088198059309244
%epsilon=1/8, Constant=12.481052357761302
%epsilon=1/16, Constant=25.085971409699642

