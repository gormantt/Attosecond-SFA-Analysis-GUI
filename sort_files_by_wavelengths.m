function [ FileName_sorted,lambda_sorted] = sort_files_by_wavelengths(PathName, FileName )
%sort_files_by_wavelengths: This function simply searches the RABBITT files
%for their wavelengths, then sorts the file names according to the
%wavelengths and returns a cell array of sorted file names by wavelength
%   Detailed explanation goes here
    lambda=zeros(length(FileName)); %declaring cell array to hold wavelengths
    for i=1:length(FileName) %searching each file for a wavelength
        fid                 =   fopen(strcat(PathName,FileName{i}),'r');

        % Header (lines starting with #)
        is_head             =   true;
        while is_head
            % Get next line
            txt             =   fgetl(fid);

            % Get the wavelength
            %   We don't want the one in eV and case
            if strmatch('#wavelength',lower(txt))
                if isempty(strmatch('#wavelength_ev',lower(txt)))
                    txt     =   strtrim(txt(strfind(txt,'=')+1:end));
                    txt     =   txt(1:strfind(txt,'m')-1);
                    lambda(i)  =   str2double(txt);   % wavelength in m     

                    txt     =   '#';
                end
            end

            % Test if still in header
            is_head         =   ( txt(1) == '#' );
        end
    end
    [lambda_sorted,sort_index]=sort(lambda);
    
    FileName_sorted=FileName(sort_index);
end

