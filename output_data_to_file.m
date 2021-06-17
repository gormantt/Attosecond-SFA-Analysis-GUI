function [] = output_data_to_file( exp_info, exp_data,det_gas,filter_name,Ip,Imin,Imax,Istep,fitmin,fitmax,filter_thickness)
%This function takes the output of SFA_main_function and the parameters
%input into SFA GUI and exports them the a file.  The parameters and such
%will appear in the header and the main data will appear in three columns:
%col 1 = integer harmonic energy (eV)
%col 2 = group delay (fs), corrected for attochirp see header for filter
%delay subtraction and/or atomic phase subtraction, currently the fit
%col 3 = group delay error (fs), currently this is the fit
%   Detailed explanation goes here
%The inputs are two cell arrays "exp_info" and "exp_data".
    for i = 1:size(exp_info,1)
        out_data=[];
        fname_for_header=char(exp_info(i,1));
        fname=fname_for_header(1:end-4);
        wavelength=cell2mat(exp_info(i,2));
        Gen_Intensity=cell2mat(exp_info(i,3));
        out_data(:,1)=transpose(cell2mat(exp_data(i,1)));
        out_data(:,2)=transpose(cell2mat(exp_data(i,2)));
        out_data(:,3)=transpose(cell2mat(exp_data(i,3)));
        fname_full = strcat('Processed_RABBITT_data\',fname,'_processed.csv');
        fid=fopen(fname_full,'w');
        fprintf(fid,'#This file has been processed through SFA_GUI \n');
        fprintf(fid,'#Source File = %s \n',fname_for_header);
        fprintf(fid,'#wavelength= %d m\n#Fitted Intensity = %d W/cm^2 \n',wavelength,Gen_Intensity);
        fprintf(fid,'#Detection Gas = %s\n#Filter=%s\n#Filter Thickness= %d nm\n',det_gas,filter_name,filter_thickness);
        fprintf(fid,'#Generation Gas Ip = %f eV\n',Ip);
        fprintf(fid,'#Intensity Fit Parameters:\n');
        fprintf(fid,'#Minimum Intensity = %fe+14 W/cm^2\n#Maximum Intensity= %fe+14 W/cm^2\n#Intensity Steps = %f\n',Imin, Imax, Istep);
        fprintf(fid,'#Fit Range Minimum = %f eV\n#Fit Range Maximum = %f eV\n', fitmin,fitmax);
        fprintf(fid,'#The ''Fit data'' option from the RABBIT_script is used for this processed data.\n');
        fprintf(fid,'#Energy (eV) |\tGroupDelay (fs) |\tError (fs)\n');
            for j=1:size(out_data,1)
                fprintf(fid,'%d\t%d\t%d\n',out_data(j,:));
            end
        
        fclose(fid);
        
    end



end

