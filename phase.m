function phi = phase(obj,t0,tr)
%phase accumulated phase for a trajectory with given ion. and rec. times
%
%   PHI = phase(T0,TR) returns the phase accumulated by a strong-field
%   approximation (SFA) trajectory with ionization and recollision times T0
%   and TR, respectively.
%
%   In the SFA, the accumulated phase is given by the integral
%      -int(T0,TR,[p + A(T)]^2/2),
%   which can be solved for ananlytically. Note that for nonzero ionization
%   potential (sfaLPshortIp and sfaLPlongIp classes) the component of the
%   phase
%      -int(T0,TR,Ip)
%   associated with the bound part of the wave function is not included. In
%   all cases, the total high-harmonic generation phase (from the stationay
%   phase equations) is simply given by
%       -int(T0,TR,[p + A(T)]^2/2 - int(T0,TR,Ip) + (ER + Ip) * TR,
%   where ER is the recollision energy (and thus ER + Ip is the harmonic
%   frequency).
%
%   See also sfaLP sfaLPshort sfaLPlong sfaLPshortIp sfaLPlongIp

%   From the definition of the lmaser field, the phase is the integral of
%      -[p - E0/omega*sin(omega*t-phi)]^2/2,
%   or equivalently
%      - (p^2 + E0^2/omega^2/2) + 2*p*E0/omega*sin(omega*t-phi) + E0^2/omega^2/2*cos(2*(omega*t-phi)),
%   which we now can integrate almost instantaneously.

% F. Mauger
%   Date            -   Revision history
%   May 24, 2016    -   Creation
%   June 9, 2016    -   Make it work for all cases 

    % Compute data
    p                   =  obj.momentum(t0,tr);
    phi                 =  -.5*(p.^2 + .5*obj.E0^2/obj.omega^2).*(tr-t0) ...
                           -obj.E0/obj.omega^2*p.*(cos(obj.omega*tr-obj.phi)-cos(obj.omega*t0-obj.phi)) ...
                           +.125*obj.E0^2/obj.omega^3*(sin(2*(obj.omega*tr-obj.phi))-sin(2*(obj.omega*t0-obj.phi)));
                              
end



