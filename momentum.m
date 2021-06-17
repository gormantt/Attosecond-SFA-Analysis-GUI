function p = momentum(obj,t0,tr)
%momentum canonical momentum for a traj. with given ion. and rec. times
%
%   P = momentum(T0,TR) returns the canonical momentum (we use the velocity
%   gauge formulation) for a strong-field approximation (SFA) trajectory
%   with ionization and recollision times T0 and TR, respectively.
%
%   From the recollision condition equation
%       p * (TR - T0) + int(T0,TR,A(T)) = 0,
%   where int is the time integration between T0 and TR, and the potential
%   vector
%       A(T) = -E0/OMEGA * sin(OMEGA*T - PHI),
%   we deduce the value for the canonical momentum P.

% F. Mauger
%   Date            -   Revision history
%   June 9, 2016    -   Creation

    % Compute data
    p                   =  -obj.E0/obj.omega^2*(cos(obj.omega*tr-obj.phi)-cos(obj.omega*t0-obj.phi))./(tr-t0);

end

