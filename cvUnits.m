classdef cvUnits
%unitConvert miscellaneous unit conversion class.
%
%Available unit convertors methods are (all static members)
%	cval = au2ev(val)
%       Energy conversion from atomic units to electron volts.
%
%   cval = au2intensity(val,ept)
%   	Conversion of the electric field amplitude E0 with ellipticity ept,
%     	expressed in atomic units, to the corresponding intensity I0,
%     	expressed in watts per squared centimeters. If unspecified, the 
%       default ellipticity corresponds to linear polarization (ept = 0).
%
%   cval = au2meter(val)
%   	Distance conversion from atomic units to meters.
%
%   cval = au2sec(val)
%   	Time conversion from atomic units to seconds.
%
%   cval = au2wavelength(val)
%     	Conversion of the electric field frequency omega, expressed in 
%       atomic units, to the corresponding wavelength lambda, expressed in
%       meters.
%
%   cval = ev2au(val)
%   	Energy conversion from electron volts to atomic units.
%
%   cval = intensity2au(val,ept)
%   	Conversion of the electric field intensity I0 with ellipticity ept,
%     	expressed in watts per squared centimeters, to the corresponding
%     	amplitude E0, expressed in atomic units. If unspecified, the 
%       default ellipticity corresponds to linear polarization (ept = 0).
%
%   cval = meter2au(val)
%   	Distance conversion from meters to atomic units.
%
%   cval = sec2au(val)
%   	Time conversion from seconds to atomic units.
%
%   cval = wavelength2au(val)
%     	Conversion of the electric field wavelength lambda, expressed in
%     	meters, to the corresponding frequency omega, expressed in atomic
%     	units.
%

% F. Mauger
%   Version 1.0.00
%   01/10/2017  Creation
    
    methods (Static=true)
%% From atomic units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cval = au2ev(val)
    % Energy conversion from atomic units to electron volts
    %   1 [a.u.] = 27.2116 [eV]
    cval                =   val*27.2116;
end
function cval = au2intensity(val,ept)
    % Conversion of the electric field amplitude E0 with ellipticity ept,
    % expressed in atomic units, to the corresponding intensity I0,
    % expressed in watts per squared centimeters
    %   I0 [W/cm²] = (1+ept^2)*(E0 [a.u.]/5.33803e-9).^2
    % If unspecified, the default ellipticity corresponds to linear
    % polarization (ept = 0).
    if nargin < 2, ept = 0; end
    cval                =   (1+ept.^2).*(val/.00000000533803).^2;
end
function cval = au2meter(val)
    % Distance conversion from atomic units to meters
    %   1 [a.u.] = 5.2917721067e-11 [m]
    cval                =   val*5.2917721067e-11;
end
function cval = au2sec(val)
    % Time conversion from atomic units to seconds
    %   1 [a.u] = 2.419e-17 [s]
    cval                =   val*2.419e-17;
end
function cval = au2wavelength(val)
    % Conversion of the electric field frequency omega, expressed in atomic
    % units, to the corresponding wavelength lambda, expressed in meters.
    %   lambda [m] = 2*pi*299792458/4.1341e16/omega [a.u.]
    cval                =   2*pi*299792458/41341000000000000./val;
end
%% To atomic units %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cval = ev2au(val)
    % Energy conversion from electron volts to atomic units
    %   1 [eV] = 1/27.2116 [a.u.]
    cval                =   val/27.2116;
end
function cval = intensity2au(val,ept)
    % Conversion of the electric field intensity I0 with ellipticity ept,
    % expressed in watts per squared centimeters, to the corresponding
    % amplitude E0, expressed in atomic units
    %   E0 [a.u.] = 5.33803e-9*sqrt(I0 [W/cm²]/(1+ept^2)).
    % If unspecified, the default ellipticity corresponds to linear
    % polarization (ept = 0).
    if nargin < 2, ept = 0; end
    cval                =   .00000000533803*sqrt(val./(1+ept^2));
end
function cval = meter2au(val)
    % Distance conversion from meters to atomic units
    %   1 [m] = 1/5.2917721067e-11 [a.u.]
    cval                =   val/5.2917721067e-11;
end
function cval = sec2au(val)
    % Time conversion from seconds to atomic units
    %   1 [s] = 1/2.419e-17 [a.u]
    cval                =   val/2.419e-17;
end
function cval = wavelength2au(val)
    % Conversion of the electric field wavelength lambda, expressed in
    % meters, to the corresponding frequency omega, expressed in atomic
    % units.
    %   omega [a.u.] = 4.1341e16/2/pi/299792458/ lambda [m]
    cval                =   2*pi*299792458/41341000000000000./val;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end

