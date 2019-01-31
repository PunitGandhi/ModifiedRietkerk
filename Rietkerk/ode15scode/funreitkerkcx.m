function Ut=funreitkerkcx(t,U,parameters)

wC=parameters.wC; %water use efficiency
Gamma=parameters.Gamma;
k1=parameters.k1;
DB0=parameters.DB0;
DW0=parameters.DW0;
kappa=parameters.kappa;
Lev=parameters.Lev;
Q=parameters.Q ;
M=parameters.M ;
VH0=parameters.VH0 ;
f=parameters.f ;

Nx=parameters.Nx ;
P0=parameters.P0;

Dx1 = parameters.Dx1;  % differentiation matrix for first derivative - not used here
Dx2 = parameters.Dx2;  % diff. matrix for second derivative

Tyear = parameters.Tyear;  % length of a year in units of equations
pseason = parameters.pseason;  % precip parameter - seasonality
Cnrm = parameters.Cnrm;        % normalization constant


tmp=reshape(U,Nx,3);

H=tmp(:,1);
W=tmp(:,2);
B=tmp(:,3);

Gw= Gamma*W./(W+k1);
In= kappa*(B + Q*f)./(B + Q);

Ht= VH0*(Dx1*H) -  In.*H + precip(t,P0,Tyear, pseason, Cnrm);  

Wt= In.*H - Lev*W - Gw.*B + DW0*(Dx2*W);

Bt= DB0*Dx2*B - M*B + wC*Gw.*B;

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

