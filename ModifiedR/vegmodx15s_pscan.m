close all
clear all
tic

global kappa Q f Sb Hd Vw Nb 
global Gamma Lev wC  phiZr Db M Kb
global Nx Dx1 Dx2 
global P Tyear alpha Cnrm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters and units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%biomass mortality
M= 0.01; %1/day
%Transpiration
Gamma=0.4000; %(cm/day)/(kg/m^2)
%water to biomass conversion
wC=0.1; %(kg/m^2)/cm
%Biomass carrying capacity
Kb=1.0; %kg/m^2
%biomass diffusion rate
Db= 0.1; %m^2/day
%soil moisture saturation
phiZr=0.45*60; %cm
%infiltration rate
Ksat= 100.0; %cm/day
%fraction of time it is actively raining
Rfrac=0.01; %0.1;%0.01;%0.001;
kappa=Ksat*Rfrac;
%infiltration contrast
f=.1;
%Biomass level for infiltration
Q= 0.1; %kg/m^2
%soil moisture evaporation rate
Lev = 0.1; %cm/day
%Surface water transport speed
Kw=4.0e5; %(m/day)/(cm)^(2/3)
Slope=0.005;
Vw=Kw*sqrt(Slope);
%Biomass dependence on surface roughness
Nb=20; %/(kg/m^2)
%exponent for surface water transport
Hd=5/3;

%Surface water height for infiltration
%HA= 0.1; %cm
%Soil moisture availability for saturation
%Sb= .05; 
Sb= 4;

%MAP=160; %mm/year
%P=MAP/365/10;
MAPvals=120:-20:40;
Pvals=MAPvals/365/10;



alpha=4;
Nacc=100000;
Cnrm=Nacc/sum(sech(alpha*cos(pi*linspace(1/Nacc,1,Nacc))).^2);
%Cnrm=12.482485478454578897;

%spatial grid, simulation parameters
TmaxD = 500*365 ; %days
dTD=5; %days between output
LxM =1000; %m
Nx = 1024;

Lx=LxM;
tmax=TmaxD;
dT=dTD;

%t=linspace(0,tmax,tmax+1);

Ntyear=round(365/dTD);
Tyear=365;



%Differentiation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx=Lx/Nx;
x = dx*(0:(Nx-1))';

dx=Lx/Nx;
x = dx*(0:(Nx-1))';


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
Sinit0=( (1- perturb/2) + perturb*intpnoise(:,2) );%.*( 1 - perturb*cos(k0x*x) );
Hinit0=( (1- perturb/2) + perturb*intpnoise(:,3) );

%Binit=ones(Nx,1);
%Sinit=ones(Nx,1);
%Hinit=ones(Nx,1);


for P=Pvals
    
close all

%%%%%%%%%%%%%%%%%%%%%%%
%Average states
%%%%%%%%%%%%%%%%%%%%%%%
%p0=T0*(80/365/10)/H0; %fix p0 so that uniform state exists
P0=P;

%bare soil
Sbs=P0/Lev;
Ibs=kappa*f*(1-Sbs)^Sb;
Hbs=P0/Ibs;
%uniform vegetation
Buv=Kb*(wC*P0-Lev*M/Gamma)/(wC*P0+M*Kb);
Suv=P0/(Lev+Gamma*Buv);
Iuv=kappa*(Buv+Q*f)/(Buv+Q) * (1-Suv)^Sb;
Huv=P0/Iuv;

    

if Buv>0
 Binit=Buv*Binit0;
 Sinit=Suv*Sinit0;
 Hinit=Huv*Hinit0;
else
 %Binit=perturb*( (1- perturb/2) + perturb*noise(:,1) ).*(1 - cos(k0x*x));
 Binit=perturb*Binit0;
 Sinit=Sbs*Sinit0;
 Hinit=Hbs*Hinit0;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0=[Hinit(:); Sinit(:); Binit(:)];
U=U0';
Ulast=U0;
tlast=0;

t=0;

ndir=0;
%outdir=['./run' int2str(ndir) '/'];
outdir=['./MAP' int2str(round(P*365*10)) ...
    'p0' int2str(perturb*100) ...
    'alpha' int2str(alpha) ...
    'kappa0' int2str(round(kappa*100)) '/'];
while exist(outdir,'dir')
    ndir=ndir+1;
    %outdir=['./run' int2str(ndir) '/'];
    outdir=[outdir(1:(end-1)) 'run' int2str(ndir) '/' ];
end
mkdir(outdir)

%store Initial Condition, Annual Averages, and Final State
fidAnnAvg= fopen([outdir 'AnnualAverages.dat'],'a+');
out=[-1, x',x',x'];
fprintf(fidAnnAvg,'%e\t',out);
fprintf(fidAnnAvg,'\n');
out=[t, U];
fprintf(fidAnnAvg,'%e\t',out);
fprintf(fidAnnAvg,'\n');


%store spatial averages, minimum, maximum for B, S H.
fidSpatAvg=fopen([outdir 'SpatialAverages.dat'],'a+');


out=[0, ...
    sum(Binit)/Nx, min(Binit), max(Binit),...
    sum(Sinit)/Nx, min(Sinit), max(Sinit),... 
    sum(Hinit)/Nx, min(Hinit), max(Hinit)];
fprintf(fidSpatAvg,'%e\t',out);
fprintf(fidSpatAvg,'\n');

Htol=1e-5;
Stol=1e-6;
Btol=1e-6;
%opts = odeset('RelTol',1e-3,'AbsTol',[Htol*ones(Nx,1);Stol*ones(Nx,1);Btol*ones(Nx,1)]);
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
        [ts, Us]=ode15s('funvegmodx',tspan,Ulast);

        %t=[t(1:end-1); ts];
        %U=[U(1:end-1,:); Us];

        Hs=Us(:,1:Nx);
        Ss=Us(:,Nx+(1:Nx));
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
        
        
        %store spatial averages, minimum, maximum for B, S H.
        out=[ts, ...
            sum(Bs,2)/Nx, min(Bs,[],2), max(Bs,[],2),...
            sum(Ss,2)/Nx, min(Ss,[],2), max(Ss,[],2),...
            sum(Hs,2)/Nx, min(Hs,[],2), max(Hs,[],2)];
        for nn=1:size(out,1)
        fprintf(fidSpatAvg,'%e\t',out(nn,:));
        fprintf(fidSpatAvg,'\n');
        end
        
        
        %print status of simulation every 10 years
        if mod(tlast/Tyear,10)<1
            display(['t=' num2str(tlast/Tyear) ' years'])
            Bdiff=sum(abs(Btmp-Ulast(2*Nx+(1:Nx))));
            Btmp=Ulast(2*Nx+(1:Nx));
            Bmin=min(Btmp);
            Bmax=max(Btmp);
            display( ['Bmin=' num2str(Bmin) ] )
            display(['Bmax='  num2str(Bmax) ] )
        end
        
        
end



fclose(fidSpatAvg);
fclose(fidAnnAvg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack output and generate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%H=Hs';%U(:,1:Nx)';
%S=Ss';%U(:,Nx+(1:Nx))';
%B=Bs';%U(:,2*Nx+(1:Nx))';
H=real(Hs');%U(:,1:Nx)';
S=real(Ss');%U(:,Nx+(1:Nx))';
B=real(Bs');%U(:,2*Nx+(1:Nx))';
t=ts;

%Save final state in simulation
out=H(1:Nx,end);
save([outdir 'Hsol.dat'],'out','-ascii')
out=S(1:Nx,end);
save([outdir 'Ssol.dat'],'out','-ascii')
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
Sanavg=AnnAvg(2:end,1+Nx+(1:Nx))';
Banavg=AnnAvg(2:end,1+2*Nx+(1:Nx))';
tanavg=AnnAvg(2:end,1);
xanavg=AnnAvg(1,1+(1:Nx))';
clear AnnAvg


%If Buv negative, flip sign for plotting to work
Bnegflag=0;
if Buv<0
  Buv=abs(Buv)
  Bnegflag=1;
end
%Spacetime plot of Annual Average B
figure1=figure();
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
box(axes1,'on');
hold(axes1,'all');
imagesc(xanavg,tanavg/Tyear,Banavg'), caxis([0,Kb/10]), colormap(flipud(summer));
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
imagesc(x,t/Tyear,B'), caxis(axB,[0,Kb/10]), colormap(axB,flipud(summer));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,round(max(t)/Tyear)])
xlim([0,Lx])
xlabel('x (m)')
ylabel('t (years)')


Hcax=Huv;
Scax=Suv;
if Hcax<0
    Hcax=Hbs;
    Scax=Hbs;
end


axS=subplot(1,3,2);
box(axS,'on');
hold(axS,'all');
imagesc(x,t/Tyear,S'), caxis(axS,[0,Scax]), colormap(axS,flipud(parula));
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
Sspavg=SpatAvg(:,5);
Hspavg=SpatAvg(:,8);
tspavg=SpatAvg(:,1);

figure('Position',[100,100,1200,600]);
pbaspect([8 1 1])
subplot(3,1,1),hold on,...
    plot(tspavg/Tyear,Bspavg,'color',[0 .5 0],'LineWidth',2),...
    plot(tspavg/Tyear,Buv*ones(size(tspavg)),'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,round(max(t)/Tyear)]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(tspavg/Tyear,Sspavg,'blue','LineWidth',1),...
    plot(tspavg/Tyear,Suv*ones(size(tspavg)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('S'),xlim([0,round(max(t)/Tyear)]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(tspavg/Tyear,Hspavg,'cyan','LineWidth',1),...
    plot(tspavg/Tyear,Huv*ones(size(tspavg)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,round(max(t)/Tyear)]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('T (years)')
%ylabel('H, S, B')
hold off
print([outdir 'WBt.eps'],'-depsc')




figure('Position',[100,100,1200,600]);
pbaspect([4 1 1])
subplot(3,1,1),hold on,...
    plot(x,Banavg(:,end),'color',[0 .5 0],'LineWidth',2),...
    plot(x,Buv*ones(size(x)),'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,Lx]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(x,Sanavg(:,end),'blue','LineWidth',1),...
    plot(x,Suv*ones(size(x)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('S'),xlim([0,Lx]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(x,Hanavg(:,end),'cyan','LineWidth',1),...
    plot(x,Huv*ones(size(x)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,Lx]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('x (m)')
%ylabel('H, S, B')
hold off
print([outdir 'HSBx.eps'],'-depsc')


%%
%%%ANIMATION
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = [outdir 'vegmod.gif'];
for nn = 1:1:length(t)
    
    % Draw plot for y = x.^n
    subplot(3,1,1),hold on,cla,...
    plot(x,B(:,nn),'color',[0 .5 0],'LineWidth',2),...
    plot(x,Buv*ones(size(x)),'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,Lx]),ylim([0,4*Buv]),pbaspect([8 1 1])
subplot(3,1,2),hold on, cla,...
    plot(x,S(:,nn),'blue','LineWidth',1),...
    plot(x,Suv*ones(size(x)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('S'),xlim([0,Lx]),ylim([0,3*Scax]),pbaspect([8 1 1])
subplot(3,1,3),hold on, cla, ...
    plot(x,H(:,nn),'cyan','LineWidth',1),...
    plot(x,Huv*ones(size(x)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,Lx]),ylim([0,25*Hcax]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('x (m)')
%ylabel('H, S, B')
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
% plot(x*X0/100,S(:,end),'b','LineWidth',2)
% plot(x*X0/100,H(:,end)*1000,'c','LineWidth',2)
% xlabel('x (m)')
% print([outdir 'HWB.png'],'-dpng')


%turn back Buv after plotting if sign was flipped
if Bnegflag
 Buv=-Buv;
end

end
