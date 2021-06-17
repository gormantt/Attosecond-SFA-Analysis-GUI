function Test
%Test Test function for sfaLP classes

%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Laser
    E0                  =   convertunits(5e13,'W/cm^2','a.u.');
    omega               =   convertunits(1700e-9,'wavelength','a.u.');
    phi                 =   0;
    
    % Target
    Ip                  =   convertunits(11.26,'eV','a.u.');
    
    % Trajectory type
    trajType            =   {'short, no Ip','long, no Ip','short','long'};
    trajCol             =   [0 1 1;1 0 1;0 0 1;1 0 0];
    
%% Test sfa predictions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization
    traj                =   sfaLP('E0',E0,'omega',omega,'phi',phi,'Ip',Ip,'type','short');
%     traj                =   sfaLP('E0',E0,'omega',omega,'phi',phi,'Ip',Ip,'type','short, no Ip');
    
    fig_time = figure('name','Ion./Rec. time');
        set(fig_time,'position',[ 365   176   659   545])
        set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
        xlabel('Ion./Rec. time (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        ylabel('Photon energy (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        
    fig_phase = figure('name','Phase');
        set(fig_phase,'position',[ 365   176   659   545])
        set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
        xlabel('Photon energy (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        ylabel('Phase (rad)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        
    fig_spectrum = figure('name','Phase');
        set(fig_spectrum,'position',[ 365   176   659   545])
        set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','log','XDir','normal','YDir','normal')
        xlabel('Photon energy (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        ylabel(' Spectrum (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        
    fig_GD = figure('name','Group delay');
        set(fig_GD,'position',[ 365   176   659   545])
        set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
        xlabel('Photon energy (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        ylabel('Group delay (rad/a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        
    fig_error = figure('name','Error');
        set(fig_error,'position',[ 365   176   659   545])
        set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
        xlabel('Photon energy (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        ylabel('Error (a.u.)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')

    % Trajectory analysis
    for k=1:length(trajType)
        % Set appropriate trajectory type
        traj.set('type',trajType{k});
        
        % Saddle point solution
        Freq          	=   linspace(1.01*traj.Ip,traj.kappa*traj.Up+traj.Ip+1,100);
        [t0,tr,p,S]     =   traj.saddlePointSolution(Freq);
        GD_app          =   real(S(2:end)-S(1:end-1))/(Freq(2)-Freq(1));
        GD              =   traj.groupDelay(Freq);
        
        % Display results
        figure(fig_time)
            hold on
                plot(real(t0),Freq,'--','Color',trajCol(k,:),'LineWidth',2)
                plot(real(tr),Freq,'-','Color',trajCol(k,:),'LineWidth',2)
            hold off
        
        figure(fig_phase)
            hold on
                plot(Freq,real(S),'-','Color',trajCol(k,:),'LineWidth',2)
            hold off
        
        figure(fig_spectrum)
            hold on
                plot(Freq,exp(-2*imag(S)),'-','Color',trajCol(k,:),'LineWidth',2)
            hold off
        
        figure(fig_GD)
            hold on
                plot(.5*(Freq(2:end)+Freq(1:end-1)),GD_app,'o','MarkerSize',10,'MarkerFaceColor','none','MarkerEdgeColor',trajCol(k,:),'LineWidth',2)
                plot(Freq,GD,'-','Color',trajCol(k,:),'LineWidth',2)
            hold off
        
        figure(fig_error)
            hold on
                plot(Freq,min(abs(traj.saddleEqIon(t0,p)-traj.Ip),abs(traj.saddleEqIon(t0,p))),':','Color',trajCol(k,:),'LineWidth',2)
                plot(Freq,abs(traj.saddleEqRec(tr,p,Freq)),'--','Color',trajCol(k,:),'LineWidth',2)
                plot(Freq,abs(traj.saddleEqTraj(t0,tr,p)),'-','Color',trajCol(k,:),'LineWidth',2)
            hold off
    end

end

