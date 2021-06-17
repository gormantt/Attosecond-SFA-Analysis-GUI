function [T,Up,te]=Trajectories(lambda,Ip,I)
% This module computes classical and semi-classical emission time
%
%     Parameters
%     ----------
%     lambda : float
%         Driving wavelength in m.
%         
%     Ip : float
%         Ionization potential in atomic units.
%         
%     I : float
%         Intensity in in 10^14 W.cm-2.
%     
%     Returns
%     -------
%     T : N x 9 matrix of floats, with N the number of odd and even
%     harmonic orders between threshold and classical cutoff.
%         T(:,1) : harmonic order.
%         T(:,2) : classical ionization time of the short trajectories.
%         T(:,3) : classical recombination time of the short trajectories.
%         T(:,4) : classical ionization time of the long trajectories.
%         T(:,5) : classical recombination time of the long trajectories.
%         T(:,6) : real part of the semi-classical ionization time of the short trajectories.
%         T(:,7) : real part of the semi-classical recombination time of the short trajectories.
%         T(:,8) : real part of the semi-classical ionization time of the long trajectories.
%         T(:,9) : real part of the semi-classical recombination time of the long trajectories.
%     Note that this module also computes the imaginary part of
%     semi-classical times but does not return them because they are
%     usually not used.
%
%     Up : float.
%         Ponderomotive potential.
%
%     te : N' x 2 matrix of floats, with N' the number of odd and even 
%     harmonic orders from threshold into the semi-classical cutoff.
%         te(:,1) : harmonic order.
%         te(:,2) : short trajectories semi-classical emission time,
%         calculated as the difference between the recombination and
%         ionization times.
%                   


%% Classical trajectories
c=299792458;
e=1.60217657e-19;
me=9.10938291e-31;
epsilon0=8.854187e-12;%
% Ip=0.01;%0.58;%in atomic units
% lambda=0.8e-6;%in meters
% I=1.4;%in 10^14W.cm-2
E0=sqrt(2*I*1e18/(epsilon0*c))/5.14e11;
% E0=0.05;%in au: 1ua =5.14e11 V./m (electric field)
omega0=(2*pi*c*24e-18./lambda);%
Up=e*(E0*5.14e11)^2*(24e-18)^2/(4*me*omega0^2);
syms t
assume(2*pi./(5*omega0)<t<2*pi./omega0)
NST=floor(4*pi./(omega0))+50;
ExcursionTime=zeros(1,NST-2);
fun=@(ti)solve(E0*(cos(omega0*t)-...
    cos(omega0*ti))/omega0^2+...
    E0/omega0*sin(omega0*ti).*(t-ti),t);
% h=waitbar(0,'Processing Classical calculations...');
for ii=3:NST
    ExcursionTime(1,ii-2)=fun((ii-1)/10);
    waitbar((ii-2)/(NST-2))
end
K=(E0*(sin(omega0*ExcursionTime)-...
    sin(omega0*(2:NST-1)/10))).^2/2/omega0^3+Ip/omega0;
StepsNb=floor(max(K)-ceil(Ip/omega0))+1;%ceil((3.17*Up/27.2+1.17*Ip)/omega0)-ceil(Ip/omega0);%floor(3.3*Up/27.2*100)+100;%
tinit=(2:NST-1)/10;
%%
en=@(ti,tr,omegaq)((E0*(sin(omega0*tr)-sin(omega0*ti)))/omega0+sqrt(2*(omegaq-Ip)));
pos=@(ti,tr)((cos(omega0*tr)-cos(omega0*ti))/omega0+sin(omega0*ti).*(tr-ti));
Ctis=zeros(1,StepsNb);Ctrs=zeros(1,StepsNb);Ctil=zeros(1,StepsNb);Ctrl=zeros(1,StepsNb);
% h3=waitbar(0,'Short Classical path...');
syms ti tr
for ii=1:StepsNb
    omegaq=(ceil(Ip/omega0)+(ii-1))*omega0;%omegaq=Ip+(ii-1)/100;%
    [Cti,Ctr]=vpasolve([en(ti,tr,omegaq),pos(ti,tr)],[ti,tr],[0,pi/omega0-10]);
    Ctis(1,ii)=Cti;
    Ctrs(1,ii)=Ctr;
    waitbar((ii-1)/(StepsNb-1))
end
%%
% h4=waitbar(0,'Long Classical path...');
for ii=1:StepsNb
    omegaq=(ceil(Ip/omega0)+(ii-1))*omega0;%((ii-1)-Ip/omega0+ceil(Ip/omega0))*omega0;%omegaq=Ip+(ii-1)/100;%
    [Cti,Ctr]=vpasolve([en(ti,tr,omegaq),pos(ti,tr)],...
        [ti,tr],[0,1.5*pi/omega0+30]);
    Ctil(1,ii)=Cti;
    Ctrl(1,ii)=Ctr;
    waitbar((ii-1)/(StepsNb-1))
end
Tinit=(ceil(Ip/omega0)+(0:size(Ctil,2)-1));
%% Quantum trajectories
p=@(Reti,Imti,Retr,Imtr)(-E0*(cos(omega0*(Retr+1i*Imtr))-...
    cos(omega0*(Reti+1i*Imti)))./(omega0^2*(Retr+1i*Imtr-(Reti+1i*Imti))));
A=@(Ret,Imt)(-E0*sin(omega0*(Ret+1i*Imt))/omega0);
fun1=@(Reti,Imti,Retr,Imtr,omegaq)(real((p(Reti,Imti,Retr,Imtr)+A(Retr,Imtr))-sqrt(2*(omegaq-Ip))));
fun2=@(Reti,Imti,Retr,Imtr,omegaq)(imag((p(Reti,Imti,Retr,Imtr)+A(Retr,Imtr))-sqrt(2*(omegaq-Ip))));
fun3=@(Reti,Imti,Retr,Imtr)(real((p(Reti,Imti,Retr,Imtr)+A(Reti,Imti))-sqrt(-2*Ip)));
fun4=@(Reti,Imti,Retr,Imtr)(imag((p(Reti,Imti,Retr,Imtr)+A(Reti,Imti))-sqrt(-2*Ip)));
syms Reti Imti Retr Imtr
ReTiqs=zeros(1,StepsNb);ImTiqs=zeros(1,StepsNb);ReTrqs=zeros(1,StepsNb);ImTrqs=zeros(1,StepsNb);
ReTiql=zeros(1,StepsNb);ImTiql=zeros(1,StepsNb);ReTrql=zeros(1,StepsNb);ImTrql=zeros(1,StepsNb);
% h1=waitbar(0,'Short quantum path...');
for ii=1:StepsNb+10
    omegaq=(ceil(Ip/omega0)+(ii-1))*omega0;%omegaq=Ip+(ii-1)/100;%
    [Imtrq,Imtiq,Retiq,Retrq]=vpasolve([fun1(Reti,Imti,Retr,Imtr,omegaq),...
        fun2(Reti,Imti,Retr,Imtr,omegaq),fun3(Reti,Imti,Retr,Imtr),fun4(Reti,Imti,Retr,Imtr)],...
        [Imtr,Imti,Reti,Retr],[0,0,0,pi/omega0-10]);
    ReTiqs(1,ii)=Retiq;
    ImTiqs(1,ii)=Imtiq;
    ReTrqs(1,ii)=Retrq;
    ImTrqs(1,ii)=Imtrq;
    waitbar((ii-1)/(StepsNb+10-1))
end
%%
% h2=waitbar(0,'Long Quantum path...');
for ii=1:StepsNb+10
    omegaq=(ceil(Ip/omega0)+(ii-1))*omega0;%omegaq=Ip+(ii-1)/100;%
    [Imtrq,Imtiq,Retiq,Retrq]=vpasolve([fun1(Reti,Imti,Retr,Imtr,omegaq),...
        fun2(Reti,Imti,Retr,Imtr,omegaq),fun3(Reti,Imti,Retr,Imtr),fun4(Reti,Imti,Retr,Imtr)],...
        [Imtr,Imti,Reti,Retr],[0,0,0,1.5*pi/omega0+30]);
    ReTiql(1,ii)=Retiq;
    ImTiql(1,ii)=Imtiq;
    ReTrql(1,ii)=Retrq;
    ImTrql(1,ii)=Imtrq;
    waitbar((ii-1)/(StepsNb+10-1))
end
close all
%%
EKQ=(ceil(Ip/omega0)+(0:StepsNb-1+10));
figure(101)
hold all
plot(K,ExcursionTime*24,'b')
plot(K,tinit*24,'b')
plot(EKQ,ReTiqs*24,'xr')
plot(EKQ,ReTrqs*24,'xr')
plot(EKQ,ReTiql*24,'xr')
plot(EKQ,ReTrql*24,'xr')
% plot(EKQ,ImTiq*24,'r--')
% plot(EKQ,ImTrq*24,'r--')
plot(Tinit,Ctil*24,'xk')
plot(Tinit,Ctis*24,'xk')
plot(Tinit,Ctrl*24,'xk')
plot(Tinit,Ctrs*24,'xk')
% plot(linspace((Ip+3.17*Up/27.2)/omega0,(Ip+3.17*Up/27.2)/omega0,2),...
%     linspace(0,2*pi/omega0*24,2),'k')
hold off
ylim([0 2*pi/omega0*24])
set(gca,'XTick',1:2:101,'XGrid','on')
%%
T(:,1)=EKQ(1:StepsNb);
T(:,2)=Ctis*24;T(:,3)=Ctrs*24;T(:,4)=Ctil*24;T(:,5)=Ctrl*24;
T(:,6)=ReTiqs(1:StepsNb)*24;T(:,7)=ReTrqs(1:StepsNb)*24;T(:,8)=ReTiql(1:StepsNb)*24;T(:,9)=ReTrql(1:StepsNb)*24;
te(:,1)=EKQ;
te(:,2)=(ReTrqs-ReTiqs)*24;