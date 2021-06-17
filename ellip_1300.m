clear;
showname;
%%
fsize=15;
%% Conversion from pixel to wavelength and Signal processing
load('wave2angle.mat');
loadname;
files=filenames;

[m n]=size(files);
NN = 2048;
wave=[-13];
wave=(wave-0.99)*2.7;
position0=1255;
halftheta=88.000;
theta0=polyval(p,-0.99*2.7);
thetam=polyval(p,wave);deltam=thetam-theta0;
theta = ((1:NN)-position0)*0.0135/800*180/pi; 
thetan=halftheta-deltam-theta;

% Generating wavelength [nm] using grating equation
%for negative value of wave, positive; for positive value of wave, negative
%22.22222: grating density
wavelength=(sind(thetan)-sind(halftheta+deltam))*(1e6/22.222222);

% Conversion from wavelength to energy
energy_p = 1239.9 ./ wavelength;
energy_calib(2:NN)= energy_p(2:NN) - energy_p(1:NN-1);
energy_calib (1) = energy_calib(2);
wave_calib(2:NN)= wavelength(2:NN) - wavelength(1:NN-1);
wave_calib(1)=wave_calib(2); wave_calib = abs(wave_calib);

% Generating signal vector
signal = zeros(m,NN);
signal_e=signal;
signal_w=signal;
exposure_time=zeros(m,1);
for ii = 1:m
    [sig back]=SifSB(files(ii,:));
    signal(ii,:) = sig;
    signal_e(ii,:)=sig./energy_calib;
    signal_w(ii,:)=sig./wave_calib;
end

%% Conversion of pressure unit
% Change the range and unit of the position according to experiments
names=load('namelist.txt');
names_l = length(names);
param=names(3:names_l);%control parameter of experiments
paramlim = [min(param)*0.99,max(param)*2];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%substract plasma line?

bg=(signal(1,:)+signal(end,:))/2;
signal=signal-ones(size(signal,1),1)*bg;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gating of individual harmonics
% plot(signal') to see spectrum vs pixel and choose the range for each HO

HH=[
400 685
690 895
900 1050
1052 1160
1162 1270
1272 1350
1352 1430
1432 1500
1502 1550
1552 1605
1607 1650
1652 1695
1697 1730
1732 1760
1762 1790
1792 1820
    ];
HHG_num=size(HH,1);
signalHH=zeros(m,HHG_num);
for ii = 1:m
    for jj = 1:HHG_num
        signalHH(ii,jj)=sum(signal(ii,HH(jj,1):HH(jj,2)));
        signalHH(ii,jj) = signalHH(ii,jj);
    end
end

%remove negative numbers in matrix "signalHH" by setting them to be eps
judge=zeros(size(signalHH));
mask=(signalHH-judge>0)-(signalHH-judge<=0)*eps;%positive number gives 1,negative gives eps;
signalHH=signalHH.*mask;
[pathstr,folder_name,ext] = fileparts(pwd);

%% Plotting
currentpath=pwd;
mkdir('Data plot');%make a folder to save pics
cd(cat(2,currentpath,'\Data plot'))

%%
%Ref. to choose HHG pixels region
sig_tot=sum(signal,2);opti=find(sig_tot==max(sig_tot));
figure,plot(signal(opti,:),'linewidth',2);
hold on
%smoothing plotting
smth_avg=5;
smth_spc=(smth_avg-1)/2;
for ii=1:NN-smth_spc
    if mod((ii+2),smth_avg)==0
        smooth(round((ii+smth_spc)/smth_avg),1)=ii;
        smooth(round((ii+smth_spc)/smth_avg),2)=sum(signal(opti,ii-smth_spc:ii+smth_spc))/smth_avg;
    end
end
plot(smooth(:,1),smooth(:,2),'r','LineWidth',2);

% Mark the range we choose to add up for each harmonic peak
xline=reshape(HH,1,numel(HH));
xline=[xline;xline];
ymax=max(signal(opti,:));
yline=[zeros(1,numel(HH));ymax*ones(1,numel(HH))];
line(xline,yline,'Color',[.8 .8 .8])

xlabel('pixel','fontsize',fsize);
ylabel('Intensity (arb. unit)','fontsize',fsize);
set(gca,'fontsize',12,'linewidth',1);
title(cat(2,num2str(param(opti)),' degree'))
ylim([-300,1.2*max(signal(opti,:))])
fig_name=cat(2,num2str(param(opti)),'degree.tif');

saveas(gcf,fig_name,'tiff');
hold off

%%
currentpath2=pwd;
mkdir('raw data energy')
cd(cat(2,currentpath2,'\raw data energy'))

for ii=1:size(signal,1)
    plot(energy_p,signal(ii,:));
    fig_name=cat(2,num2str(param(ii)),' degree.tif');
    title(fig_name)
    saveas(gcf,fig_name,'tiff');
end

cd(currentpath2)
%%
currentpath2=pwd;
mkdir('raw data')
cd(cat(2,currentpath2,'\raw data'))

for ii=1:size(signal,1)
    plot(signal(ii,:));
    fig_name=cat(2,num2str(param(ii)),' degree.tif');
    title(fig_name)
    saveas(gcf,fig_name,'tiff');
end

cd(currentpath2)

%% Looking for angle00
angle=names(3:end);
fit_y=[];waist=[];angle0=[];peak=[];fitting_n=[];fit_y_n=[];
angle_int=linspace(min(angle),max(angle),5*length(angle))
figure
for ii=1:size(signalHH,2)
fit_y=signalHH(:,ii);
%x is fitting parameter array; X0 is initial parameter array
X0=[max(fit_y),angle(find(fit_y==max(fit_y))),10];
x=lsqnonlin(@(x,X,Y) x(1)*exp(-(X-x(2)).^2/x(3)^2) - Y,X0,[],[],[],angle,fit_y);
angle0=[angle0,x(2)];
fitting=x(1)*exp(-(angle_int-x(2)).^2/x(3).^2);

plot(angle,fit_y,'*',angle_int,fitting,'k-')
    fig_name=cat(2,num2str(ii),' (HO).tif');
    title(fig_name)
    saveas(gcf,fig_name,'tiff');
end
angle0
HO_upper=7
angle00=mean(angle0(1:HO_upper))
%%
cd(currentpath)%back to original folder
