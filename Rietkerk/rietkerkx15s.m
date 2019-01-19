close all
clear all
tic

% global Dh Dw  
% global a f p g m v
% global Nx
% global Dx1 Dx2 
% global p0 Tyear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters and units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=10; %conversion by plant
gmax=0.05; 
k1=5; %mm
d=0.05; %/day
DB0=0.1;%m^2/day
fr=0.005; %raining 2 percent of the time 
alpha=0.2; %/day
k2=5; %g/m^2
W0r=0.2;
rw=0.2;%/day
DW0=0.1;%m^2/day
DH0=100;%m^2/day

MAP=450;% mm/year
R=MAP/365; %mm/day


LxM =500; %m
Nx = 512;

parameters.c = c;
parameters.gmax = gmax;
parameters.k1 = k1;
parameters.d = d;
parameters.DB0 = DB0;
parameters.fr = fr;
parameters.alpha = alpha;
parameters.k2 = k2;
parameters.W0r = W0r;
parameters.rw = rw;
parameters.DW0 = DW0;
parameters.DH0 = DH0;
parameters.R = R;
parameters.LxM = LxM;
parameters.Nx = Nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nondimensial parameters
%%%%%%%%%%%%%%%%%%%%R=30/3*24; %mm/day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k1=k1/fr;
H0=k1; %mm
W0=k1; %mm
B0=k2; %g/m^2
T0=1/(c*gmax); %days
X0=sqrt(DB0*T0); %m

Lx=LxM/X0;

Dw=DW0/DB0;
Dh=DH0/DB0;
a=alpha*T0;
f=W0r;
m=d*T0;
g=k2/(k1*c);
v=rw*T0;
p0=R*T0/k1;


parameters.Dw = Dw;
parameters.Dh = Dh;
parameters.a = a;
parameters.f = f;
parameters.m = m;
parameters.g = g;
parameters.v = v;
parameters.p0 = p0;



%Differentiation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx=Lx/Nx;
x = dx*(0:(Nx-1))';


%finite difference differentiation matrices with PBC
%Dx1=(diag(ones(Nx-1,1),1)-diag(ones(Nx-1,1),-1) )/(2*dx); 
%Dx1(1,Nx)=-1/(2*dx); Dx1(Nx,1)=1/(2*dx);%PBC

e = ones(Nx,1);

%first derivative with upwinding
Dx1 = spdiags([-e e],[0 1],Nx,Nx); 
%PBC
Dx1(Nx,1)=1;

%centered second derivative
Dx2=spdiags([e -2*e e],-1:1,Nx,Nx); 
%PBC
Dx2(Nx,1)=1; Dx2(1,Nx)=1;


Dx1=Dx1./dx;
Dx2=Dx2./dx.^2;

parameters.Dx1 = Dx1;
parameters.Dx2 = Dx2;

%initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(462)
noise=rand(Nx,2);
perturb=.1;



Binit=( (1- perturb/2) + perturb*noise(:,1) );%.*( 1 - perturb*cos(k0x*x) );
Winit=( (1- perturb/2) + perturb*noise(:,2) );%.*( 1 - perturb*cos(k0x*x) );
Hinit=0*noise(:,2); ( (1- perturb/2) + perturb*noise(:,2) );


Wuv= m/(1-m);
Buv= (p0 - v*Wuv) / (g*m);
Iuv=a*fr*(Buv+f)/(Buv+1);
Huv=p0/Iuv;

Wbs=p0/v;
Hbs=p0/(a*f*fr);

if Buv>0
 Binit=Buv*Binit;
 Winit=Wuv*Winit;
 Hinit=Huv*Hinit;
else
 %Binit=perturb*( (1- perturb/2) + perturb*noise(:,1) ).*(1 - cos(k0x*x));
 Binit=perturb*Binit;
 Winit=Wbs*Winit;
 Hinit=Hbs*Hinit;
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0=[Hinit(:); Winit(:); Binit(:)];
U=U0';
Ulast=U0;
tlast=0;

t=0;

ndir=0;
outdir=['./run' int2str(ndir) '/'];
%outdir=['./MAP' int2str(MAP) 'p' int2str(perturb*100) '/'];
while exist(outdir,'dir')
    ndir=ndir+1;
    outdir=['./run' int2str(ndir) '/'];
    %outdir=[outdir(1:10) 'run' int2str(ndir) '/' ];
end
mkdir(outdir)



Tday=1/T0;
Thour=Tday/24;
Tmin=Thour/60;
Tyear=Tday*365;

parameters.Tyear = Tyear;

tspan=linspace(0,500*Tyear,101);

        
        [tr, Ur]=ode15s(@(t,U) funreitkerkcx(t,U,parameters),tspan,Ulast);

        Hr=Ur(:,1:Nx);
        Wr=Ur(:,Nx+(1:Nx));
        Br=Ur(:,2*Nx+(1:Nx));
        
       Ulast=Ur( end , 1:(3*Nx) );
        %Ulast(1:Nx)=Ulast(1:Nx)+p;
        tlast=tr(end);
            

        
        
        
        
    
figure('Position',[100,100,1200,600]);
pbaspect([8 1 1])
subplot(3,1,1),hold on,...
    plot(tr*T0,B0*sum(Br,2)/Nx,'color',[0 .5 0],'LineWidth',2),...
    ylabel('B (g/m^2)'),xlim([0,tr(end)*T0]),ylim([0,1.5]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(tr*T0,sum(Wr,2)/Nx*W0,'blue','LineWidth',1),...
    ylabel('W (mm)'),xlim([0,tr(end)*T0]),ylim([0,30]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(tr*T0,H0*sum(Hr,2)/Nx,'cyan','LineWidth',1),...
    ylabel('H (mm)'),xlim([0,tr(end)*T0]),ylim([0,30]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('T (days)')
%ylabel('H, S, B')
hold off
print(['WBtr.eps'],'-depsc')


figure('Position',[100,100,1200,400],'DefaultAxesFontSize',18);
set(gca,'FontSize',18) 
pbaspect([16 1 1])
subplot(3,1,1),hold on,...
    plot(tl*T0,B0*sum(Bl,2)/Nx,'color',[0 .5 0],'LineWidth',2),...
    ylabel('B (g/m^2)','FontSize',18),xlim([0,tl(end)*T0]),ylim([0,1.5]),pbaspect([16 1 1])
    fill([0,.5,.5,0,],[0 0 1.5 1.5],[0,0,1],'FaceAlpha',0.1)
subplot(3,1,2),hold on,...
    plot(tl*T0,sum(Wl,2)/Nx*W0,'blue','LineWidth',1),...
    ylabel('W (mm)','FontSize',18),xlim([0,tl(end)*T0]),ylim([0,35]),pbaspect([16 1 1])
    fill([0,.5,.5,0,],[0 0 35 35],[0,0,1],'FaceAlpha',0.1)
subplot(3,1,3),hold on,...
    plot(tl*T0,H0*sum(Hl,2)/Nx,'cyan','LineWidth',1),...
    ylabel('H (mm)','FontSize',18),xlim([0,tl(end)*T0]),ylim([0,20]),pbaspect([16 1 1])
    fill([0,.5,.5,0,],[0 0 20 20],[0,0,1],'FaceAlpha',0.1)
%ylim([0,1])
xlabel('T (days)','FontSize',18)
%ylabel('H, S, B')
hold off
print(['WBtl.png'],'-dpng')
