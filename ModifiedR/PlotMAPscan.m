close all
clear all
tic
global v d p Bk Ha Sb
global eps sig rho mu  q f
global gamma Ia
global Nx Lx
global Dx1 Dx2 
global Tyear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters and units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%biomass mortality
m= .005; %1/day
%Transpiration
Gamma=2000; %cm^3/kg/day
%water to biomass conversion
w=0.1; %1/cm
%Biomass carrying capacity
Kb=2e-4; %kg/cm^2

%soil moisture saturation
phiZr=30; %cm
%infiltration rate
K= 100; %cm/day
%infiltration contrast
f=.1;
%Biomass level for infiltration
Q= 10^-6; %kg/cm^2
%Surface water height for infiltration
HA= 1.0; %cm
%Soil moisture availability for saturation
Sb= .05; 
%soil moisture evaporation rate
r = 0.1; %cm/day

MAP=300; %mm/year
%precipitation rate during rainfall event
%R=4; %cm/day

%roughness parameter
Nb=5e4; %80000;

Kw=3e8;
Db=2;
Slope=0.005;
%X0=1000; %cm Spatial lengthscale

%spatial grid, simulation parameters
TmaxD = 200*365 ; %days
dTD=5; %days between output
LxM =500; %m
Nx = 512;

%t=linspace(0,tmax,tmax+1);

Ntyear=round(365/dTD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nondimensial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0=HA;
B0=Kb;
T0=H0/K;
Slope0=sum(Slope)/length(Slope);
X0=Kw*T0*sqrt(Slope0);
%X0=10000;
Lx=LxM*100/X0;

eps=w*Gamma*Kb*T0;
mu = m*T0/eps;

Bk = Kb/B0;

q   =Q/B0;
Ha  =HA/H0;
sig =r/(Gamma*B0);
gamma= 1/(w*phiZr);
Ia= K*T0/phiZr;

rho=B0*Nb;

p0=T0*(MAP/365/10)/H0;
p=p0; %use constant avg rainfall

v=sqrt(Slope)*Kw*T0/X0; %v includes factor of \sqrt(Slope/Slope0) 
d=T0*Db/X0^2/eps;
%d=T0*Db/X0^2;

tmax=TmaxD/T0;
dT=dTD/T0;
%%%%%%%%%%%%%%%%%
%Rain parameters%
%%%%%%%%%%%%%%%%%
%Tseason=100/T0;
%Nrain=5;
Tyear=365/T0;


%Tperiod=Tseason/Nrain;

%dimensionless diffusion
%d=1; 
%dimensionless advection
%v=5000;

%Slope
%Sg=0.005; 
%water flow constant
%Kw= 1; 
%Biomass diffusion rate
%Db=1;

%%%%%%%%
%set lengthscale
%X0=sqrt(Db*phiZr/K); %cm to make d=1

%Mean Annual Precipitation
%MAP=p*H0*(Nrain+1); %cm/year

%p0=MAP/K/365;


%%%%%%%%%%%%%%%%%%%%%%%
%Average states
%%%%%%%%%%%%%%%%%%%%%%%
%bare soil
Sbs=p0*Ia/eps/sig/gamma;
Ibs=f*(1-Sbs)/(1+Sb-Sbs);
Hbs=Ha*p0/(Ibs-p0);
%uniform vegetation
Buv=Bk*(p0*Ia-eps*gamma*mu*sig)/(p0*Ia+eps*gamma*mu*Bk);
Suv=(p0*Ia+eps*gamma*mu*Bk)/(eps*gamma)/(Bk+sig);
Iuv=(Buv+q*f)/(Buv+q) * (1-Suv)/(1-Suv+Sb);
Huv=Ha*p0/(Iuv-p0);
 





%Differentiation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx=Lx/Nx;
x = dx*(0:(Nx-1))';


%finite difference differentiation matrices with PBC
%Dx1=(diag(ones(Nx-1,1),1)-diag(ones(Nx-1,1),-1) )/(2*dx); 
%Dx1(1,Nx)=-1/(2*dx); Dx1(Nx,1)=1/(2*dx);%PBC


%centered second derivative with PBC
Dx2=(diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1) - 2*diag(ones(Nx,1),0) )/(dx^2);
Dx2(1,Nx)=1/(dx^2); Dx2(Nx,1)=1/(dx^2);%PBC


%finite difference differentiation matrices with upwinding
Dx1=(diag(ones(Nx-1,1),1)-diag(ones(Nx,1),0) )/(dx);
Dx1(Nx,1)=1/(dx);%PBC






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAPvals=[20, 30, 40, 50, 60,70,80,90,100,105,110,120,130];
p0vals=T0*(MAPvals/365/10)/H0;
BUVvals=Bk*(p0vals*Ia-eps*gamma*mu*sig)./(p0vals*Ia+eps*gamma*mu*Bk);
figure1=figure;
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
box(axes1,'on');
hold(axes1,'all');
plot(MAPvals,BUVvals*B0*1e4*1000,'Color',[0 .5 0],'LineStyle','--','LineWidth',2)
%plot stable UV states
MAPUVs=[105,110,120];
BUVs= [2.065217e-01, 2.200855e-01, 2.458678e-01];
plot(MAPUVs,BUVs*B0*1e4,'s','MarkerSize',12,'MarkerFace','k')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAPvals=[20, 30, 40, 50, 60,70,80,90,100,101,105,110,120,130];
p0vals=T0*(MAPvals/365/10)/H0;
BUVvals=Bk*(p0vals*Ia-eps*gamma*mu*sig)./(p0vals*Ia+eps*gamma*mu*Bk);
figure1=figure;
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[2 1 1],'FontSize',16);
box(axes1,'on');
hold(axes1,'all');

%predicted Turing-Hopf
MAPTH=101;
nunst=MAPvals<MAPTH;
%plot(101*[1 1],[0,1],'b:')
%plot predicted UV states
plot(MAPvals(nunst),BUVvals(nunst)*B0*1e4,'Color',[0 .5 0],'LineStyle','--','LineWidth',2)
plot(MAPvals(~nunst),BUVvals(~nunst)*B0*1e4,'Color',[0 .5 0],'LineStyle','-','LineWidth',2)

%plot max and avg for patterns
patcol=[0.5 0.5 0.5];
MAPP=[30, 40, 50, 60,70,80,90,100,105];
BPavg=[5.474925e-02,	7.612436e-02, 9.752556e-02, 1.189549e-01, 1.404133e-01,1.617640e-01,1.830787e-01, 2.040169e-01, 2.096552e-01];
BPmax=[	3.077625e-01, 3.083431e-01, 3.107003e-01, 3.154685e-01, 3.175172e-01,3.161657e-01,3.114882e-01, 3.025714e-01, 2.688336e-01];
plot(MAPP,BPavg*B0*1e4,'o','MarkerSize',12,'MarkerEdge',patcol,'MarkerFace',patcol)
%plot(MAPP,BPmax*B0*1e4,'^','MarkerSize',12,'MarkerEdge',patcol,'MarkerFace',patcol)


for nn=1:length(MAPP)
    plot(MAPP(nn)*[1 1],[0,BPmax(nn)]*B0*1e4,'Color',patcol)
end


%plot stable UV states
MAPUVs=[105,110,120];
BUVs= [2.065217e-01, 2.200855e-01, 2.458678e-01];
plot(MAPUVs,BUVs*B0*1e4,'s','MarkerSize',12,'MarkerFace','k','MarkerEdge','k')

xlim([20 130])
ylim([0 0.7])
xlabel('MAP (mm/year)')
ylabel('B (kg/m^2)')



pulsewidth=zeros(5,1);
pulseheight=zeros(5,1);
pulseheightW=zeros(5,1);
pulsep0=zeros(5,1);
dirs=dir('MAP*p10');
nc=0;
for nd=1:length(dirs)
    B=load(['./' dirs(nd).name '/Bsol.dat']);
    tmp=load(['./' dirs(nd).name '/sol.mat']);
    tmpBuv=abs(tmp.Buv);
    pltflg=0;
    if min(B)< tmpBuv/5;
        pltflg=1;
    end
    
    nn=1;
    pbcflag=0;
    nstp=0;
    
    if pltflg
        %check that pulse goes across boundary
        stph=0;
        while B(nn)>tmpBuv && nn<Nx;
            nn=nn+1;
            pbcflag=1;
            stph=max(nstp,(B(nn)));
        end
        nstp=nn;
        
        
        while nn<Nx
            tmpw=0;
            tmph=0;
            if B(nn)< tmpBuv
                nn=nn+1;
            else
                nc=nc+1;
                while B(nn)>tmpBuv && nn<Nx
                    nn=nn+1;
                    tmpw=tmpw+1;
                    tmph=max(tmph,B(nn));
          
                end
                if nn==Nx
                    tmpw=tmpw+nstp;
                    tmph=max(tmph,stph);
                end
                pulsewidth(nc)=tmpw*dx;
                pulseheight(nc)=tmph;
                pulsep0(nc)=tmp.MAP;
            end
        end
    end
end
figure1=figure;
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[2 1 1],'FontSize',16);
box(axes1,'on');
hold(axes1,'all');
gscatter(pulsewidth(pulsep0~=105)*X0/100,pulseheight(pulsep0~=105)*B0*1e4,pulsep0(pulsep0~=105))
%gscatter(pulsewidth*X0/100,pulseheight*B0*1e4,pulsep0)
%xlim([40 80])
%ylim([0.579 0.64])
xlabel('pulse width (m)')
ylabel('B_{max} (kg/m^2)')



xxx

















%initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(462)
noise=rand(Nx,2);
%B0 = (a/2/m+real(sqrt((a/2/m)^2-1)))*(0.9995+0.001*noise(:,1));
%W0 = m*(a/2/m-real(sqrt((a/2/m)^2-1)))*(0.9995+0.001*noise(:,2));
%m=mu;
%a0=.1;
k0x=5*(2*pi/Lx);
%B0 = (a0/2/m+real(sqrt((a0/2/m)^2-1)))*ones(Nx,1).*(1+cos(k0x*x));
%W0 = m*(a0/2/m-real(sqrt((a0/2/m)^2-1)))*ones(Nx,1).*(1-cos(k0x*x));
%H0= 0*ones(Nx,1);
%filter

Binit=Buv*(0.95+0.1*noise(:,1)).*(1-.1*cos(k0x*x));
Sinit=Suv*(0.95+0.1*noise(:,2)).*(1-.1*cos(k0x*x));
Hinit=Huv*(0.95+0.5*noise(:,2));

%Binit=Bk*ones(Nx,1);
%Sinit=zeros(Nx,1);
%Hinit=zeros(Nx,1);

%Bhat = fft(B0);
%B0=real(ifft(Bhat));
%What=fft(W0);
%W0=real(ifft(W0));







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0=[Hinit(:); Sinit(:); Binit(:)];
U=U0';
Ulast=U0;
tlast=0;

t=0;

ndir=0;
outdir=['./run' int2str(ndir) '/'];
while exist(outdir,'dir')
    ndir=ndir+1;
    outdir=['./run' int2str(ndir) '/'];
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
opts = odeset('RelTol',1e-3,'AbsTol',[Htol*ones(Nx,1);Stol*ones(Nx,1);Btol*ones(Nx,1)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run Simulation year by year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while tlast<tmax
     
        tspan=tlast+(0:dT:Tyear);
        
        [ts, Us]=ode15s('funvegmodx',tspan,Ulast,opts);
        %t=[t(1:end-1); ts];
        %U=[U(1:end-1,:); Us];

        Hs=Us(:,1:Nx);
        Ss=Us(:,Nx+(1:Nx));
        Bs=Us(:,2*Nx+(1:Nx));
        
        
        
        Ulast=Us( end , 1:(3*Nx) );
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
            display( ['Bmin=' num2str(min(Ulast(2*Nx+(1:Nx)))) ] )
            display(['Bmax='  num2str(max(Ulast(2*Nx+(1:Nx)))) ] )
        end
        
        
end



fclose(fidSpatAvg);
fclose(fidAnnAvg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack output and generate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=Hs';%U(:,1:Nx)';
S=Ss';%U(:,Nx+(1:Nx))';
B=Bs';%U(:,2*Nx+(1:Nx))';
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


%Spacetime plot of Annual Average B
figure1=figure();
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',16);
box(axes1,'on');
hold(axes1,'all');
imagesc(xanavg*X0/100,tanavg/Tyear,Banavg'), caxis([0,2.5*Buv]), colormap(flipud(summer));
set(gca,'YDir','normal')
ylim([0,tmax/Tyear])
xlim([0,Lx*X0/100])
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
imagesc(x*X0/100,t/Tyear,B'), caxis(axB,[0,2.5*Buv]), colormap(axB,flipud(summer));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,tmax/Tyear])
xlim([0,Lx*X0/100])
xlabel('x (m)')
ylabel('t (years)')


axS=subplot(1,3,2);
box(axS,'on');
hold(axS,'all');
imagesc(x*X0/100,t/Tyear,S'), caxis(axS,[0,Suv]), colormap(axS,flipud(parula));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,tmax/Tyear])
xlim([0,Lx*X0/100])
xlabel('x (m)')
ylabel('t (years)')

axH=subplot(1,3,3);
box(axH,'on');
hold(axH,'all');
imagesc(x*X0/100,t/Tyear,H'), caxis(axH,[0,Huv]), colormap(axH,flipud(parula));
set(gca,'YDir','normal')
ylim([min(t)/Tyear,tmax/Tyear])
xlim([0,Lx*X0/100])
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
    plot(tspavg/Tyear,B0*Bspavg*10^4,'color',[0 .5 0],'LineWidth',2),...
    plot(tspavg/Tyear,B0*Buv*ones(size(tspavg))*10^4,'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,tmax/Tyear]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(tspavg/Tyear,Sspavg,'blue','LineWidth',1),...
    plot(tspavg/Tyear,Suv*ones(size(tspavg)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('S'),xlim([0,tmax/Tyear]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(tspavg/Tyear,H0*Hspavg,'cyan','LineWidth',1),...
    plot(tspavg/Tyear,H0*Huv*ones(size(tspavg)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,tmax/Tyear]),pbaspect([8 1 1])

%ylim([0,1])
xlabel('T (years)')
%ylabel('H, S, B')
hold off
print([outdir 'WBt.eps'],'-depsc')




figure('Position',[100,100,1200,600]);
pbaspect([4 1 1])
subplot(3,1,1),hold on,...
    plot(x*X0/100,B0*Banavg(:,end)*10^4,'color',[0 .5 0],'LineWidth',2),...
    plot(x*X0/100,B0*Buv*ones(size(x))*10^4,'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,Lx*X0/100]),pbaspect([8 1 1])
subplot(3,1,2),hold on,...
    plot(x*X0/100,Sanavg(:,end),'blue','LineWidth',1),...
    plot(x*X0/100,Suv*ones(size(x)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('S'),xlim([0,Lx*X0/100]),pbaspect([8 1 1])
subplot(3,1,3),hold on,...
    plot(x*X0/100,H0*Hanavg(:,end),'cyan','LineWidth',1),...
    plot(x*X0/100,H0*Huv*ones(size(x)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,Lx*X0/100]),pbaspect([8 1 1])

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
    plot(x*X0/100,B0*B(:,nn)*10^4,'color',[0 .5 0],'LineWidth',2),...
    plot(x*X0/100,B0*Buv*ones(size(x))*10^4,'color',[0 .5 0],'LineWidth',2,'LineStyle',':'),...  
    ylabel('B (kg/m^2)'),xlim([0,Lx*X0/100]),ylim([0,4*Buv*B0*10^4]),pbaspect([8 1 1])
subplot(3,1,2),hold on, cla,...
    plot(x*X0/100,S(:,nn),'blue','LineWidth',1),...
    plot(x*X0/100,Suv*ones(size(x)),'blue','LineWidth',1,'LineStyle',':'),...    
    ylabel('S'),xlim([0,Lx*X0/100]),ylim([0,3*Suv]),pbaspect([8 1 1])
subplot(3,1,3),hold on, cla, ...
    plot(x*X0/100,H0*H(:,nn),'cyan','LineWidth',1),...
    plot(x*X0/100,H0*Huv*ones(size(x)),'cyan','LineWidth',1,'LineStyle',':'),...    
    ylabel('H (cm)'),xlim([0,Lx*X0/100]),ylim([0,25*Huv*H0]),pbaspect([8 1 1])

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
