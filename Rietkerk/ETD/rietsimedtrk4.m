function [Tout, Uout] = rietsimedtrk4(Uinit,tspan,parameters)


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

Lx=parameters.LxM;
Nx=parameters.Nx ;
P0=parameters.P0;
dT=parameters.dTD;

%Nout=parameters.Nout;

Tyear = parameters.Tyear;  % length of a year in units of equations
pseason = parameters.pseason;  % precip parameter - seasonality
Cnrm = parameters.Cnrm;        % normalization constant


tmp=reshape(Uinit,Nx,3);

H=tmp(:,1);
W=tmp(:,2);
B=tmp(:,3);



x = Lx*(1:Nx)'/Nx;

%wavenunmber space
kx = 2*pi*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx; % wave numbers


%Define Linear part of equation
rietLH = 1i*VH0*kx;
rietLW = -DW0*kx.^2 - Lev;
rietLB = -DB0*kx.^2 - M; 



[Eh, E2h, Qh,f1h,f2h,f3h] = etdrk4coeffs(rietLH,dT);
[Ew, E2w, Qw,f1w,f2w,f3w] = etdrk4coeffs(rietLW,dT);
[Eb, E2b, Qb,f1b,f2b,f3b] = etdrk4coeffs(rietLB,dT);




t=tspan(1);
Uout(1,:)= Uinit;
Tout(1)=t;

nmax = round(diff(tspan)/dT);

%filter 1/3 of wavenumbers for (approximate) quadratic nonlinearity
xfilter=(floor(Nx/3)+1):(ceil(2*Nx/3)+1);

[Hhat,What,Bhat]=takeFFT(H,W,B,xfilter);


for nn = 2:(nmax+1)
    t = t+dT;
    
    H0hat=Hhat;
    W0hat=What;
    B0hat=Bhat;
    
    %%step 1
    %%%%%%%%   
   [NLh0,NLw0,NLb0]=rietNLhat(H,W,B,t,parameters,xfilter); 
   H1hat= E2h.*H0hat+Qh.*NLh0;
   W1hat= E2w.*W0hat+Qw.*NLw0;
   B1hat= E2b.*B0hat+Qb.*NLb0;

   %step 2
   %%%%%%%
   [Ht,Wt,Bt]=takeIFFT(H1hat,W1hat,B1hat);
   [NLh1,NLw1,NLb1]=rietNLhat(Ht,Wt,Bt,t,parameters,xfilter); 
   H2hat= E2h.*H0hat+Qh.*NLh1;
   W2hat= E2w.*W0hat+Qw.*NLw1;
   B2hat= E2b.*B0hat+Qb.*NLb1;
   
   %step 3
   %%%%%%%%
   [Ht,Wt,Bt]=takeIFFT(H2hat,W2hat,B2hat);
   [NLh2,NLw2,NLb2]=rietNLhat(Ht,Wt,Bt,t,parameters,xfilter); 
   H3hat= E2h.*H1hat+Qh.*(2*NLh2-NLh0);
   W3hat= E2w.*W1hat+Qw.*(2*NLw2-NLw0);
   B3hat= E2b.*B1hat+Qb.*(2*NLb2-NLb0);
   
   %final output
   %%%%%%%%
   [Ht,Wt,Bt]=takeIFFT(H3hat,W3hat,B3hat);
   [NLh3,NLw3,NLb3]=rietNLhat(Ht,Wt,Bt,t,parameters,xfilter); 
   Hhat= Eh.*H0hat + f1h.*NLh0 + 2*f2h.*(NLh1+NLh2) + f3h.*NLh3;
   What= Ew.*W0hat + f1w.*NLw0 + 2*f2w.*(NLw1+NLw2) + f3w.*NLw3;
   Bhat= Eb.*B0hat + f1b.*NLb0 + 2*f2b.*(NLb1+NLb2) + f3b.*NLb3;
   
   [H,W,B]=takeIFFT(Hhat,What,Bhat);
   Uout(nn,:)=[H;W;B]';
   Tout(nn)=t;
end


    
function [Hhat,What,Bhat]=takeFFT(H,W,B,xfilter)

Hhat=fft(H);
Hhat(xfilter)=0;

What=fft(W);
What(xfilter)=0;

Bhat=fft(B);
Bhat(xfilter)=0;

function [H,W,B]=takeIFFT(Hhat,What,Bhat)

H=real(ifft(Hhat));
W=real(ifft(What));
B=real(ifft(Bhat));

function [NLh,NLw,NLb]=rietNLhat(H,W,B,t,parameters,xfilter)

wC=parameters.wC;
kappa=parameters.kappa;
Q=parameters.Q;
f=parameters.f;
Gamma=parameters.Gamma;
k1=parameters.k1;
Tyear=parameters.Tyear;
pseason=parameters.pseason;
Cnrm=parameters.Cnrm;
P0=parameters.P0;


Infl=kappa*(B+Q*f)./(B+Q);
Grwth=Gamma*W./(W+k1);

NLh=fft(precip(t,P0,Tyear, pseason, Cnrm)-Infl.*H);
NLh(xfilter)=0;

NLw=fft(Infl.*H-Grwth.*B);
NLw(xfilter)=0;

NLb=fft(wC*Grwth.*B);
NLb(xfilter)=0;



%% Precipitation function

function pt=precip(t,p0,Tyear, pseason, Cnrm)
%pt=p0  *  sech(  cos(pi*t/Tyear) / epsilon ).^2  * Constant
 pt=p0  *  sech(pseason*sin(pi*t/Tyear)           ).^2  *  Cnrm;
 
%these values may only be accurate to 3 or 4 decimals.
%epsilon=1, Constant=1.494515826022178
%epsilon=1/2, Constant=2.806475267441896
%epsilon=1/4, Constant=6.088198059309244
%epsilon=1/8, Constant=12.481052357761302
%epsilon=1/16, Constant=25.085971409699642