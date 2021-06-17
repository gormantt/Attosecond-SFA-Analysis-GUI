function [exported_info,exported_data] =SFA_main_function_v2(Imin,Imax,Istep,Ip,fitmin,fitmax,filter_name,filter_thickness,det_gas)

    %Title: SFA_main_function
    %Author: Tim Gorman
    %This function loads RABBITT data created by Dietrich's RABBIT script and it
    %subtracts a filter delay, an atomic delay from detection or neither
    %then it subtracts the group delay calculated by SFA.
    %The last version of the script "SFA_script" used before creating this
    %function was "SFA_script_v3".  This function also outputs the
    %subtracted data and the parameters used in it's subtraction.
    %v2


    [FileName,PathName,FilterIndex] = uigetfile(fullfile(pwd,'RABBITT_data','*.csv'), 'Pick a RABBITT file!','MultiSelect','on');
    addpath(genpath(pwd));
    %close all;

    if length(FileName)>1 %checks whether or not files were selected
        fig_compare=figure('name','Plot All with Errorbars');
        set(fig_compare,'position',[ 900  176   780   545]); %figure with error bars
        set(gca,'Box','on','LineWidth',0.5,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal');
        xlabel('Photon Energy (eV)','FontSize',20,'FontName','Times New Roman','FontAngle','normal');
        ylabel('Group delay (fs)','FontSize',20,'FontName','Times New Roman','FontAngle','normal');
        title('All with Attosub','FontSize',20,'FontName','Times New Roman','FontAngle','normal');
        fig_compare2=figure('name','Plot All without Errorbars'); %figure without error bars
        set(fig_compare2,'position',[ 900  176   780   545]);
        set(gca,'Box','on','LineWidth',0.5,'FontSize',20,'XScale','linear','YScale','linear','XDir','normal','YDir','normal');
        xlabel('Photon Energy (eV)','FontSize',20,'FontName','Times New Roman','FontAngle','normal');
        ylabel('Group delay (fs)','FontSize',20,'FontName','Times New Roman','FontAngle','normal');
        title('All with Attosub','FontSize',20,'FontName','Times New Roman','FontAngle','normal');

        %Defualt parameters that were used in script version before making
        %this a function.
        %Imin=0.25;% 10^14 W/cm^2
        %Imax=1.0;% 10^14 W/cm^r
        %Istep=40; %steps between Imin and Imax
        %Ip=12.61; %Ionization potential of Generation Gas.
        %fitmin=25; %eV
        %fitmax=50; %eV
        %filter_name='none';
        %filter_thickness=200;
        %det_gas='none';
        
        exported_info={};%information about data to export from function
        exported_data={};%data to export from function
        I0_m_list={};
        
        fitrange=[fitmin, fitmax];%range over which SFA is fit
        legendinfo={};
        filt_phase  = filter_phase(filter_name,filter_thickness); %loading filter phase for subtraction

        [at_delay, det_Ip]=load_atomic_delay(det_gas); %delay (as) vs. electron kinetic energy (eV) and the detection atomic ip (eV)
        at_delay_mcc=at_delay_minus_cc(at_delay); % extrapolated and interpolated delay (as) vs. electron kinetic energy (eV) with 800 cc subtracted


        %Check to see if one file or multiple files

        if iscell(FileName)==true %there are multiple files

            [FileName_sorted,lambda_sorted]=sort_files_by_wavelengths(PathName,FileName); %sort file names according to wavelength in RABBITT file
            exported_info(:,1)=FileName_sorted(:,1);
            exported_info(:,2)=num2cell(lambda_sorted(:,1));
            col  = winter(length(FileName));
            for i=1:length(FileName)

                    loopfile=strcat(PathName,FileName_sorted{i});
                    [E_GD,GD_sub,GD_exp_err,lambda,I0_m]=intensityFit_v4(loopfile,Imin,Imax,Istep, Ip, fitrange,filt_phase,at_delay_mcc,det_Ip); %this is where the fitting is done

                    figure(fig_compare)
                    hold on
                    errorbar(E_GD,GD_sub,GD_exp_err,'x','Color',col(i,:),'LineWidth',2)
                    hold off
                    drawnow
                    figure(fig_compare2)
                    hold on
                    plot(E_GD,GD_sub,'x','Color',col(i,:),'LineWidth',2)
                    hold off
                    legendinfo{i}=[num2str(lambda*10^9)];
                    I0_m_list(i,1)=num2cell(I0_m);
                    exported_data(i,1)=mat2cell(E_GD,1);%energy axes to be exported by function
                    exported_data(i,2)=mat2cell(GD_sub,1);%Group Delay data to be exported by function
                    exported_data(i,3)=mat2cell(GD_exp_err,1); %Group delay error to be exported by function
            end
            exported_info(:,3)=I0_m_list;
        elseif iscell(FileName)==false && ~isempty(FileName) %there is only one file
            loopfile=strcat(PathName,FileName);
            exported_info(:,1)=cellstr(FileName);
            [E_GD,GD_sub,GD_exp_err,lambda,I0_m]=intensityFit_v4(loopfile,Imin,Imax,Istep, Ip, fitrange,filt_phase,at_delay_mcc,det_Ip); %this is where the fitting is done

            figure(fig_compare)
            hold on
            errorbar(E_GD,GD_sub,GD_exp_err,'x','LineWidth',2)
            hold off
            drawnow
            figure(fig_compare2)
            hold on
            plot(E_GD,GD_sub,'x','LineWidth',2)
            hold off
            legendinfo{1}=[num2str(lambda*10^9)];
            exported_info(:,2)=num2cell(lambda);
            exported_info(:,3)=num2cell(I0_m);
            exported_data(1,1)=mat2cell(E_GD,1);%energy axes to be exported by function
            exported_data(1,2)=mat2cell(GD_sub,1);%Group Delay data to be exported by function
            exported_data(1,3)=mat2cell(GD_exp_err,1);%Group delay error to be exported by function
        end
        %xlabel('Photon Energy (eV)','FontSize',20,'FontName','Times New Roman','FontAngle','normal');
        figure(fig_compare);
        legend(legendinfo,'Location','NorthEastOutside');

        figure(fig_compare2);
        legend(legendinfo,'Location','NorthEastOutside');
    else %no files chosen
        sprintf('You did not choose a file!')
    end
end
 
 
