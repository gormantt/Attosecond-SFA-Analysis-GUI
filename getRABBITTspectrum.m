function [E,Pos,RAB] = getRABBITTspectrum(fname,units)
%getRABBITTspectrum return the RABBITT spectrum from input file
%
%   [E,Pos,RAB] = getRABBITTspectrum(FileName) reads and returns the
%   RABBITT trace spectrum from the FileName file. Output vectors E and Pos
%   are the energies of the harmonics and position of the motor,
%   respectively. RAB is a length(Pos)-by-length(E) matrix of the RABBITT
%   spectrum (linear are for  a constant position and rows for a constant
%   energy).
%
%   [...] = getRABBITTspectrum(FileName,Units) specifies the Units to be
%   used for Energy values
%       'phys'/'eV' uses electron volts,
%       'a.u.'      uses atomic units.
%   Default Units is 'phys'. Additional unit conversions can be performed
%   with the cvUnits class.
%
%   See also cvUnits

% F. Mauger
%   Version 1.0.00
%   01/11/2017  Creation

%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Units
    if nargin < 2, units = []; end
    
    if isempty(units)
        units           =   'phys';
    elseif ~any(strcmpi(units,{'a.u.','au','atomicunits','atomic_units','phys','physical','real','ev','electronvolts','electron_volts'}))
        warning('getRABBITTspectrum:Units',['Unknown set of units ' units '. Use default phys instead'])
        units           =   'phys';
    end
    
    % Output display
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%                                         %%%')
    disp('%%%          RABBITT spectrum data          %%%')
    disp('%%%                                         %%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Read data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open file
    fid                 =   fopen(fname,'r');
    disp(fname)
    disp(' ')
    
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
                
                disp(['   Wavelength        ' deblank(txt) ' m'])
                
                txt     =   '#';
            end
        elseif strmatch('#detection gas',lower(txt))
            txt         =   strtrim(txt(strfind(txt,'=')+1:end));
                
                disp(['   Detection gas     ' deblank(txt)])
                
                txt     =   '#';
        elseif strmatch('#ip',lower(txt))
            txt         =   strtrim(txt(strfind(txt,'=')+1:end));
            txt         =   txt(1:strfind(txt,'e')-1);
                
                disp(['   Ip                ' deblank(txt) ' eV'])
                
                txt     =   '#';
        end
        
        % Test if still in header
        is_head         =   ( txt(1) == '#' );
    end
    
    % Actual data
    D                   =   [];
    while ischar(txt)
        % Convert data
        D               =   [D;str2num(txt)];
        
        % Get next set of data
        txt             =   fgetl(fid);
    end
    
    % Close file
    fclose(fid);
    disp(' ')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
%% Extract relevant data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E                   =   D(1,2:end);
    Pos                 =   D(2:end,1).';
    RAB                 =   D(2:end,2:end);
    
    if any(strcmpi(units,{'a.u.','au','atomicunits','atomic_units'}))
        E               =   cvUnits.ev2au(E);
    end
    
end

