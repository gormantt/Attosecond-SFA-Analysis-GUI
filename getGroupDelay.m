function [E,GD1,GD2,GD3,GD4] = getGroupDelay(fname,units)
%getGroupDelay return the group delay from an input RABBITT trace file
%
%   [E,GD1,GD2,GD3,GD4] = getGroupDelay(FileName) reads and returns the
%   group energy E and group delays GD1-4 from the RABBITT file FileName.
%   The group delays are computed from various RABBITT fits
%       GD1     FFT (From Scan)
%       GD2     FFT (From Int Scan)
%       GD3     Fit (From Scan)
%       GD4     Fit (From Int Scan)
%
%   [...] = getGroupDelay(FileName,Units) specifies the Units to be used in
%   output results. Available units are
%       'phys'  returns the energy in electron volts and group delay in
%               femtoseconds,
%       'a.u.'  returns both the energy and group delay in atomic units.
%   Default system of units is 'phys'. Additional unit conversions can be
%   performed with the cvUnits class.
%
%   See also cvUnits

% F. Mauger
%   Version 1.0.00
%   01/07/2017  Creation
%   01/11/2017  Add unit selection option

%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Units
    if nargin < 2, units = []; end
    
    if isempty(units)
        units           =   'phys';
    elseif ~any(strcmpi(units,{'a.u.','au','atomicunits','atomic_units','phys','physical','real'}))
        warning('getGroupDelay:Units',['Unknown set of units ' units '. Use default phys instead'])
        units           =   'phys';
    end

%% Read data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open file
    fid                 =   fopen(fname,'r');
    
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
                lambda  =   str2double(txt);
                
                txt     =   '#';
            end
        end
        
        % Test if still in header
        is_head         =   ( txt(1) == '#' );
    end
    
    % Laser frequency
    omega               =   cvUnits.wavelength2au(lambda);
    if any(strcmpi(units,{'phys','physical','real'}))
        omega           =   cvUnits.sec2au(omega)*1e-15;
    end
    
    % Actual data
    GD                  =   [];
    while ischar(txt)
        % Convert to data
        G               = str2num(txt);
        
        % Update results >>> we only keep even harmonics
        % Data file organization (rows)         -|
        %   1   Harmonic Order (int)             |
        %   2   Harmonic Order (actual)          |
        %   3   Harmonic Energy (eV)             |  Energy and Harmonic Order
        %   4   Window Start (index)             |
        %   5   Window End (index)              _|
        %   6   RABBITT Phase (Unwrapped) (rad) -|  FFT (From Scan)
        %   7   RABBITT Phase Error (rad)       _|
        %   8   RABBITT Phase (Unwrapped) (rad) -|  FFT (From Int Scan)
        %   9   RABBITT Phase Error (rad)       _|
        %   10  RABBITT Phase (Unwrapped) (rad) -|  Fit (From Scan)
        %   11  RABBITT Phase Error (rad)       _|
        %   12  RABBITT Phase (Unwrapped) (rad) -|  Fit (From Int Scan)
        %   13  RABBITT Phase Error (rad)       _|
        if mod(G(1),2)==0
            GD          =   [GD; G([3 6 8 10 12])];                         %#ok<AGROW>
        end
        
        % Get next set of data
        txt             =   fgetl(fid);
    end
    
    % Close file
    fclose(fid);

%% Compute group delays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GD                  =   GD.';
    E                   =   GD(1,:);
    GD1                 =   GD(2,:)*.5/omega;
    GD2                 =   GD(3,:)*.5/omega;
    GD3                 =   GD(4,:)*.5/omega;
    GD4                 =   GD(5,:)*.5/omega;
    
    if any(strcmpi(units,{'a.u.','au','atomicunits','atomic_units'}))
        E               =   cvUnits.ev2au(E);
    end
    
end