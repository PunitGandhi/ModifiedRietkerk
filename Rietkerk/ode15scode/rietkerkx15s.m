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
wC=0.1; %(kg/m^2)/cm  conversion by plant
Gamma=5; %(cm/day)/(kg/m^2) 
k1=0.5; %cm
DB0=0.1;%m^2/day
kappa=0.2; %0.2; 100 %/day
Lev=0.2; %/day
M=0.25;%/day
Q= 0.005; %kg/m^2
DW0=0.1;%m^2/day
VH0= 10; %10; 3e4;%m/day
f=0.1;

MAP=400;% mm/year
P0=MAP/365/10; %cm/day


LxM =500; %m
Nx = 512;

Tmax=500*365; %days
dTD=5; %days

parameters.wC = wC; %water use efficiency
parameters.Gamma = Gamma;
parameters.k1 = k1;
parameters.DB0 = DB0;
parameters.DW0 = DW0;
parameters.kappa = kappa;
parameters.Lev = Lev;
parameters.Q = Q;
parameters.M = M;
parameters.VH0 = VH0;
parameters.f = f;
parameters.LxM = LxM;
parameters.Nx = Nx;
parameters.P0=P0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nondimensial parameters
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0=1;
W0=1;
B0=1;
T0=1;
X0=1;
%H0=k1; %mm
%W0=k1; %mm ***currently can't change
%B0=Q; %g/m^2 ***currently can't change
%T0=1/(kappa); %days
%X0=VH0*T0; %m

Tday=1/T0;
Thour=Tday/24;
Tmin=Thour/60;
Tyear=Tday*365;


Lx=LxM/X0;
tmax= Tmax/Tday;
dT=dTD/T0;

eps=wC*Gamma*T0;

Dw=T0*DW0/X0^2/eps;
Db=T0*DB0/X0^2/eps;

mu=M*T0/eps;
gamma=B0*T0*Gamma/W0/eps;
sig=Lev*T0/eps;

alpha=kappa*T0;
zeta=H0/W0;
nu=VH0*T0/X0;

p0=P0*T0/H0;


pseason=0;
Nacc=100000;
Cnrm=Nacc/sum(sech(pseason*cos(pi*linspace(1/Nacc,1,Nacc))).^2);

% parameters.Dw = Dw;
% parameters.Db = Db;
% parameters.eps= eps;
% parameters.f = f;
% parameters.mu = mu;
% parameters.gamma = gamma;
% parameters.sig = sig;
% parameters.p0 = p0;
% 
% parameters.alpha=alpha;
% parameters.zeta=zeta;
% parameters.nu=nu;


parameters.Tyear = Tyear;
parameters.pseason = pseason;
parameters.Cnrm = Cnrm;



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
%Limit high frequency of noise
rng(462)
Nxn=1024; 
Lxn=5000;
noise=rand(Nxn,3);
xnoise = (Lxn/Nxn)*(0:(Nxn-1))';
intpnoise=interp1(xnoise,noise,x,'pchip');
%B0 = (a/2/m+real(sqrt((a/2/m)^2-1)))*(0.9995+0.001*noise(:,1));
%W0 = m*(a/2/m-real(sqrt((a/2/m)^2-1)))*(0.9995+0.001*noise(:,2));
%m=mu;
%a0=.1;
%k0x=4*(2*pi/Lx);
%B0 = (a0/2/m+real(sqrt((a0/2/m)^2-1)))*ones(Nx,1).*(1+cos(k0x*x));
%W0 = m*(a0/2/m-real(sqrt((a0/2/m)^2-1)))*ones(Nx,1).*(1-cos(k0x*x));
%H0= 0*ones(Nx,1);
%filter

perturb=.01;

Binit0=( (1- perturb/2) + perturb*intpnoise(:,1) );%.*( 1 - perturb*cos(k0x*x) );
Winit0=( (1- perturb/2) + perturb*intpnoise(:,2) );%.*( 1 - perturb*cos(k0x*x) );
Hinit0=( (1- perturb/2) + perturb*intpnoise(:,3) );



W0uv= k1*M/(wC*Gamma-M);
B0uv= wC*(P0 - Lev*W0uv) / M;
I0uv=kappa*(B0uv+Q*f)/(B0uv+Q);
H0uv=P0/I0uv;

W0bs=P0/Lev;
H0bs=P0/(kappa*f);

if B0uv>0
 Binit=B0uv*Binit0/B0;
 Winit=W0uv*Winit0/W0;
 Hinit=H0uv*Hinit0/H0;
else
 %Binit=perturb*( (1- perturb/2) + perturb*noise(:,1) ).*(1 - cos(k0x*x));
 Binit=perturb*Binit0;
 Winit=W0bs*Winit0/W0;
 Hinit=H0bs*Hinit0/H0;
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
%outdir=['./run' int2str(ndir) '/'];
outdir=['./MAP' int2str(round(P0*365*10)) ...
    'alpha' int2str(pseason) ...
    'kappa' int2str(round(kappa*100)) ...
    'V' int2str(VH0) '/'];
while exist(outdir,'dir')
    ndir=ndir+1;
    %outdir=['./run' int2str(ndir) '/'];
    if ndir>1
    outdir=[outdir(1:(end-5)) 'run' int2str(ndir) '/' ];
    else
     outdir=[outdir(1:(end-1)) 'run' int2str(ndir) '/' ];
    end
end
mkdir(outdir)
display([outdir])

%store Initial Condition, Annual Averages, and Final State
fidAnnAvg= fopen([outdir 'AnnualAverages.dat'],'a+');
out=[-1, x',x',x'];
fprintf(fidAnnAvg,'%e\t',out);
fprintf(fidAnnAvg,'\n');
out=[t, U];
fprintf(fidAnnAvg,'%e\t',out);
fprintf(fidAnnAvg,'\n');


%store spatial averages, minimum, maximum for B, W H.
fidSpatAvg=fopen([outdir 'SpatialAverages.dat'],'a+');


out=[0, ...
    sum(Binit)/Nx, min(Binit), max(Binit),...
    sum(Winit)/Nx, min(Winit), max(Winit),... 
    sum(Hinit)/Nx, min(Hinit), max(Hinit)];
fprintf(fidSpatAvg,'%e\t',out);
fprintf(fidSpatAvg,'\n');

Htol=1e-5;
Wtol=1e-6;
Btol=1e-6;
%opts = odeset('RelTol',1e-3,'AbsTol',[Htol*ones(Nx,1);Wtol*ones(Nx,1);Btol*ones(Nx,1)]);
opts=odeset('Stats','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run Simulation year by year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bdiff=1;
Bmax=max(Binit);
eps0=1e-10;
Btmp=Binit';
while tlast<tmax && Bdiff>eps0 && Bmax>eps0
        %display(num2str(tlast));
        tspan=tlast+(0:dT:Tyear);
        
        %[ts, Us]=ode15s('funvegmodx',tspan,Ulast,opts);
        %[ts, Us]=ode15s('funreitkerkcx',tspan,Ulast);
        [ts, Us]=ode15s(@(t,U) funreitkerkcx(t,U,parameters),tspan,Ulast);
        %toc
        %t=[t(1:end-1); ts];
        %U=[U(1:end-1,:); Us];

        Hs=Us(:,1:Nx);
        Ws=Us(:,Nx+(1:Nx));
        Bs=Us(:,2*Nx+(1:Nx));
        
        
        Ulast= Us( end , 1:(3*Nx) );
        %Ulast=real(Us( end , 1:(3*Nx) ));
        %Ulast(1:Nx)=Ulast(1:Nx)+p;
        tlast=ts(end);
        
        
        Uannavg=sum(Us,1)/length(ts);
        %store Annual Averages
        out=[tlast, Uannavg];
        fprintf(fidAnnAvg,'%e\t',out);
        fprintf(fidAnnAvg,'\n');
        
        
        %store spatial averages, minimum, maximum for B, W H.
        out=[ts, ...
            sum(Bs,2)/Nx, min(Bs,[],2), max(Bs,[],2),...
            sum(Ws,2)/Nx, min(Ws,[],2), max(Ws,[],2),...
            sum(Hs,2)/Nx, min(Hs,[],2), max(Hs,[],2)];
        for nn=1:size(out,1)
        fprintf(fidSpatAvg,'%e\t',out(nn,:));
        fprintf(fidSpatAvg,'\n');
        end
        
        
        %print status of simulation every 10 years
        if mod(tlast/Tyear,10)<1
            display(['t=' num2str(tlast/Tyear) ' years'])
            Bdiff=sum(abs(Btmp-Ulast(2*Nx+(1:Nx))))/Nx/10;
            Btmp=Ulast(2*Nx+(1:Nx));
            Bmin=min(Btmp);
            Bmax=max(Btmp);
            display( ['Bmax=' num2str(Bmax) ] )
            display(['Bamplitude='  num2str(Bmax-Bmin) ] )
            display(['Bdiff='  num2str(Bdiff) ] )
            toc
        end
        
        
end



fclose(fidSpatAvg);
fclose(fidAnnAvg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack output and generate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%H=Hs';%U(:,1:Nx)';
%W=Ws';%U(:,Nx+(1:Nx))';
%B=Bs';%U(:,2*Nx+(1:Nx))';
H=real(Hs');%U(:,1:Nx)';
W=real(Ws');%U(:,Nx+(1:Nx))';
B=real(Bs');%U(:,2*Nx+(1:Nx))';
t=ts;

%Save final state in simulation
out=H(1:Nx,end);
save([outdir 'Hsol.dat'],'out','-ascii')
out=W(1:Nx,end);
save([outdir 'Wsol.dat'],'out','-ascii')
out=B(1:Nx,end);
save([outdir 'Bsol.dat'],'out','-ascii')


%save workspace including final year of data
save([outdir 'sol.mat'])


toc

%%
%load annual averages for spacetime plot
AnnAvg=load([outdir 'AnnualAverages.dat']);
 %skip 1st row with x values and 1st column with t values
 %note that first row of data will be initial condition
Hanavg=AnnAvg(2:end,1+(1:Nx))';
Wanavg=AnnAvg(2:end,1+Nx+(1:Nx))';
Banavg=AnnAvg(2:end,1+2*Nx+(1:Nx))';
tanavg=AnnAvg(2:end,1);
xanavg=AnnAvg(1,1+(1:Nx))';
clear AnnAvg


%If B0uv negative, flip sign for plotting to work
Bnegflag=0;
if B0uv<0
  B0uv=abs(B0uv)
  Bnegflag=1;
end
%Spacetime plot of Annual Average B
figure1=figure();
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
box(axes1,'on');
hold(axes1,'all');
imagesc(xanavg,tanavg/Tyear,Banavg'), caxis([0,B0uv]), colormap(flipud(summer));
set(gca,'YDir','normal')
ylim([0,round(tanavg(end)/Tyear)])
xlim([0,Lx])
xlabel('x (m)')
ylabel('t (years)')
print([outdir 'Bspacetime.png'],'-dpng')


%%
%Spacetime plot of B,S,H over last year of simulation
figure1=figure('Position',[100,100,1200,300]);
axB=subplot(1,3,1);
%axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
box(axB,'on');
hold(axB,'all');
imagesc(x,t/Tyear,B'), caxis(axB,[0,B0uv]), colormap(axB,flipud(summer));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,round(max(t)/Tyear)])
xlim([0,Lx])
xlabel('x (m)')
ylabel('t (years)')


Hcax=H0uv;
Wcax=W0uv;
if Hcax<0
    Hcax=Hbs;
    Wcax=Hbs;
end


axW=subplot(1,3,2);
box(axW,'on');
hold(axW,'all');
imagesc(x,t/Tyear,W'), caxis(axW,[0,Wcax]), colormap(axW,flipud(parula));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,round(max(t)/Tyear)])
xlim([0,Lx])
xlabel('x (m)')
ylabel('t (years)')

axH=subplot(1,3,3);
box(axH,'on');
hold(axH,'all');
imagesc(x,t/Tyear,H'), caxis(axH,[0,Hcax]), colormap(axH,flipud(parula));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,round(max(t)/Tyear)])
xlim([0,Lx])
xlabel('x (m)')
ylabel('t (years)')

print([outdir 'spacetime1.png'],'-dpng')


%%

%Load spatial averages
SpatAvg=load([outdir 'SpatialAverages.dat']);
Bspavg=SpatAvg(:,2);
Wspavg=SpatAvg(:,5);
Hspavg=SpatAvg(:,8);
tspavg=SpatAvg(:,1);

figure('Position',[100,100,1200,600]);
pbaspect([8 1 1])
subplot(3,1,1),hold on,...
    plot(tspavg/Tyear,Bspavg,'color',[0 .5 0],'LineWidth',2),...
    plot(tspavg/Tyear,B0uv*ones(size(tspavg)),'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,round(max(t)/Tyear)]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(tspavg/Tyear,Wspavg,'blue','LineWidth',1),...
    plot(tspavg/Tyear,W0uv*ones(size(tspavg)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('W'),xlim([0,round(max(t)/Tyear)]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(tspavg/Tyear,Hspavg,'cyan','LineWidth',1),...
    plot(tspavg/Tyear,H0uv*ones(size(tspavg)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,round(max(t)/Tyear)]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('T (years)')
%ylabel('H, W, B')
hold off
print([outdir 'WBt.eps'],'-depsc')




figure('Position',[100,100,1200,600]);
pbaspect([4 1 1])
subplot(3,1,1),hold on,...
    plot(x,Banavg(:,end),'color',[0 .5 0],'LineWidth',2),...
    plot(x,B0uv*ones(size(x)),'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,Lx]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(x,Wanavg(:,end),'blue','LineWidth',1),...
    plot(x,W0uv*ones(size(x)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('W'),xlim([0,Lx]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(x,Hanavg(:,end),'cyan','LineWidth',1),...
    plot(x,H0uv*ones(size(x)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,Lx]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('x (m)')
%ylabel('H, W, B')
hold off
print([outdir 'HWBx.eps'],'-depsc')


%%
%%%ANIMATION
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = [outdir 'vegmod.gif'];
for nn = 1:1:length(t)
    
    % Draw plot for y = x.^n
    subplot(3,1,1),hold on,cla,...
    plot(x,B(:,nn),'color',[0 .5 0],'LineWidth',2),...
    plot(x,B0uv*ones(size(x)),'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,Lx]),ylim([0,4*B0uv]),pbaspect([8 1 1])
subplot(3,1,2),hold on, cla,...
    plot(x,W(:,nn),'blue','LineWidth',1),...
    plot(x,W0uv*ones(size(x)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('W'),xlim([0,Lx]),ylim([0,3*Wcax]),pbaspect([8 1 1])
subplot(3,1,3),hold on, cla, ...
    plot(x,H(:,nn),'cyan','LineWidth',1),...
    plot(x,H0uv*ones(size(x)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,Lx]),ylim([0,25*Hcax]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('x (m)')
%ylabel('H, W, B')
hold off

    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if nn == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.1); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1); 
      end
end

  

% figure1=figure();
% axes1 = axes('Parent',figure1,'FontSize',16);
% box(axes1,'on');
% hold(axes1,'all');
% plot(x*X0/100,B(:,end),'color',[0,.5,0],'LineWidth',2)
% plot(x*X0/100,W(:,end),'b','LineWidth',2)
% plot(x*X0/100,H(:,end)*1000,'c','LineWidth',2)
% xlabel('x (m)')
% print([outdir 'HWB.png'],'-dpng')


%turn back B0uv after plotting if sign was flipped
if Bnegflag
 B0uv=-B0uv;
end


