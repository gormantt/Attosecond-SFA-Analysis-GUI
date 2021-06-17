function [E_GD,GD_sub,GD_exp_err,lambda,I0_m]=intensityFit_v4(fname_in,Imin,Imax,Istep,Ip_in,fit_range,filt_phase,at_delay_mcc,det_Ip)
%intensityFit Fit the effective intensity for group delay measurement
%v2: adds functionality to externally interact with the function
%v3: adds functionality to add filter and atomic phase subtraction before fitting for
%intensity
%v4: exports the all of the intensities that minimized the fits to a single
%matrix.  These will be used in the GUI if the files are to be exported
    
%% Parameters and Read Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input file
    fname               =   fname_in;%'2016-08-16_19-48-46_2d_avg.csv';
    fGroupDelay         =   [fname];
    %Load filter data and calculate phase
    %filt_phase          = filter_phase(filter_name,filter_thickness);
    % Target parameters
    Ip                  =   cvUnits.ev2au(Ip_in); % Ip_in is in eV
    hbar=6.582119514*10^-16; % eV*s/2pi
    tol=0.001;  %tolerance in eV for aluminum filter range calculations
    
    [E_GD,~,~,~,~,~, GD_exp, GD_exp_err, ~, ~,lambda]       =   getGroupDelay_with_error(fGroupDelay,'phys');%Currently using integer value for harmonic energy and extracting the fit
    
    fprintf('wavelength=%d\n',lambda)
    
    
    %wavelength in eV
    lambda_eV=1240/(lambda*10^9);
    
    %setting up x-axes involved in calculating GD due to filter
    odd_harm_E=lambda_eV:2*lambda_eV:1000*lambda_eV;
    even_harm_E=2*lambda_eV:2*lambda_eV:1000*lambda_eV;

    %Subtracting Filter GD
    %Interpolate filter phase to odd harmonic energies from the 1st to the 1000th harmonic and calculate filter
    %delay if the filter_phase value is not empty
    if isempty(filt_phase)~=1
        inter_filt_phase=interp1(filt_phase(:,1),filt_phase(:,2),odd_harm_E,'spline');

        filter_delay=diff(inter_filt_phase)./diff(odd_harm_E)*hbar*10^15; %delay in fs relative to vacuum
        even_harm_E_delay=even_harm_E(1:length(filter_delay));
        
        %Find edges with tolerance listed above
        even_start=even_harm_E_delay(abs(E_GD(1)-even_harm_E_delay)<=tol);
        even_end=even_harm_E_delay(abs(E_GD(end)-even_harm_E_delay)<=tol);

        filt_range=(even_harm_E_delay >= even_start) & (even_harm_E_delay<=even_end);
        
        fprintf('E_GD length =%d\n',length(E_GD));
        fprintf('Even_harm=%d\n',length(even_harm_E_delay(filt_range)));
        
        GD_exp=GD_exp-filter_delay(filt_range);       
    end
    
    %Subtracting atomic delay
    if isempty(at_delay_mcc)~=1 && isempty(det_Ip)
        at_cc_delay=cc_delay(lambda);
        at_delay_mcc(:,2)=at_delay_mcc(:,2)+at_cc_delay(:,2); %adding back in CC delay for this wavelength
        at_delay_mcc(:,1)=at_delay_mcc(:,1)+det_Ip;%changing from electron KE to photon energy
        at_delay_int=interp1(at_delay_mcc(:,1),at_delay_mcc(:,2),E_GD,'spline'); %interpolate to the experimental x-axis
        GD_exp=GD_exp-at_delay_int;
    end
    

    % Laser parameters
    omega               =   cvUnits.wavelength2au(lambda); % lambda is in m
    I0                  =   linspace(Imin,Imax,Istep)*1e14;%input intensity is in units of 10^14 W/cm^2
    E0                  =   cvUnits.intensity2au(I0);
    
    % Data analysis
    GD_fit              =   fit_range;
    GD_shift            =   GD_fit;
    
    % Display parameters
    col                 =   jet(length(I0));
    
    fig_GD = figure('name','Group delay');
        set(fig_GD,'position',[ 150   176   659   545])
        set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
        hold on
            % Experimental group delay
            errorbar(E_GD,GD_exp,GD_exp_err,'k-','LineWidth',2)
            
            % Fit region
            ylim(ylim())
            plot(GD_fit(1)*[1 1],ylim(),'k--','LineWidth',1)
            plot(GD_fit(2)*[1 1],ylim(),'k--','LineWidth',1)
        hold off
        xlabel('Photon energy (eV)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        ylabel('Group delay (fs)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
        
%% Intensity scan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization
    E                   =   cvUnits.ev2au(E_GD);
    fs                  =   cvUnits.au2sec(1)*1e15;
    traj                =   sfaLP('E0',E0(1),'omega',omega,'phi',0,'Ip',Ip,'type','short','refFreq',E);
    
    Err                 =   NaN(size(I0));
    ind_shift           =   (E_GD >= GD_shift(1)) & (E_GD <= GD_shift(2));
    ind_fit             =   (E_GD >= GD_fit(1)) & (E_GD <= GD_fit(2));

    
    % Intensity scan
    for k=1:length(I0)
        % Update trajectory
        traj.set('E0',E0(k)); %E0(k) is the E-field intensity
        
        % Group delay computation and error
        GD              =   traj.groupDelay(E)*fs;
        %GD              =   GD + mean(GD_exp(ind_shift)-GD(ind_shift));
        shift_guess     =   mean(GD(ind_shift)-GD_exp(ind_shift));
        
        GD        =   GD - simpleshiftfit(GD(ind_shift),GD_exp(ind_shift),[],GD_exp_err(ind_shift),shift_guess);
        
        Err(k)          =   sqrt(sum(((GD(ind_fit)-GD_exp(ind_fit)).^2)./GD_exp_err(ind_shift).^2));
        
        % Display on curve
        figure(fig_GD)
            hold on
                plot(E_GD,GD,'--','Color',col(k,:),'LineWidth',2)
                %errorbar(E_GD,GD_exp,GD_exp_err,'k-','LineWidth',2)
            hold off
            drawnow
        
    end
        
    [~,ind]             =   min(Err);
    I0_m                =   I0(ind);
    
    I                   =   griddedInterpolant(I0,Err,'pchip','pchip');
    I0_m                =   fminsearch(@(i) I(i),I0_m); % in units of 10^14 W/cm^2
    sprintf('I0_m = %e',I0_m);
    I0_m_string = num2str(I0_m,'%e');
    title_string = strcat('Minimum Intensity = ', I0_m_string);
    figure(fig_GD)
        title(title_string)
    
    % Error curve
%     fig_Err = figure('name','Error');
%         set(fig_Err,'position',[ 365   176   659   545])
%         set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
%         hold on
%             plot(I0,Err,'r-','LineWidth',2)
%             plot(I0_m*[1 1],ylim,'k--','LineWidth',2)
%         hold off
%         xlabel('Intensity (W/cm^{2})','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
%         ylabel('Group delay error (fs)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
%         title(['Best fit I0 ' num2str(I0_m,3) ' W/cm^{-2}'],'FontSize',20,'FontName','Times New Roman','FontAngle','normal')

    %Subtracted Curve
    E0_min = cvUnits.intensity2au(I0_m);
    traj.set('E0',E0_min);
    GD_min=traj.groupDelay(E)*fs;
    shift_guess_min     =   mean(GD_min(ind_shift)-GD_exp(ind_shift));
    GD_sub = (GD_exp+simpleshiftfit(GD_min(ind_shift),GD_exp(ind_shift),[],GD_exp_err(ind_shift),shift_guess_min))-GD_min;
    
%     fig_Sub=figure('name','subtract')
%     set(fig_Sub,'position',[ 365   176   659   545])
%     set(gca,'Box','on','LineWidth',2,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal')
%         hold on
%         errorbar(E_GD,GD_sub,GD_exp_err);
%         hold off
%     xlabel('Photon Energy (eV)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
%     ylabel('Group delay (fs)','FontSize',20,'FontName','Times New Roman','FontAngle','normal')
%     title(['Attochirp Subtracted, I0 = ' num2str(I0_m,3) ' W/cm^{-2}'],'FontSize',20,'FontName','Times New Roman','FontAngle','normal')

end

