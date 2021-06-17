classdef (Abstract) sfaLPsolver
%Interface class for linear polarization strong-field approximation solver.
%This class is part of the sfaLP suite and is not meant to be used on its
%own by end users.
    
%% Interface properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Abstract,GetAccess=public,SetAccess=private)
        type                                % Trajectory type
    end
    
    properties (GetAccess=public,SetAccess=protected)
        t0_c                                % Ionization time for maximum return energy
        t0_m                                % Ionization time for minimum return energy
        tr_c                                % Recollision time for maximum return energy
        tr_m                                % Recollision time for minimum return energy
        kappa                               % Maximum return energy (Up)
    end
    
%% Nonliner solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=protected,Static=true,Sealed=true)
function X = NLsolver(varargin)
    % Nonlinear system solver to be actually used
    X                   =   systemSolve.fsolve(varargin{:});%    loc_fsolve(varargin{:});
end
    end
    
%% Interface methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Abstract)
        [t0,tr,p,S]     =   saddlePointSolution(obj,sfaParam,freq)
    end
    
%% Initialization methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=protected)
function obj = sfaLPsolver(sfaParam,trajIndex)
    % Constructor
    
if nargin > 1
    % Cut-off trajectory
    % The cut-off trajectory is obtained from the stationary phase Eqs.
    % [(3-5) in the sfaLP documentation] ignoring the ionization potential
    % in the ionization condition [Ip=0 in (3)]. From there we deduce the
    % canonical momentum and recollision energy
    %   p = -A(t0)      and         Er = (p+A(tr))^2/2,
    % while satisfying the recollision condition
    %   p*(t_r-t_0) + int(A(s),s=t_0,t_r) = 0.
        % Function handles for cut-off trajectory determination
        traj            =   @(t0,t) sfaParam.saddleEqTraj(t0,t,-sfaParam.potVector(t0));
        tr              =   @(t0) fzero(@(tr) traj(t0,tr), (4.399-sfaParam.phi)/sfaParam.omega);
        
        % Appropriate initial condition
        if (trajIndex == 1) || (trajIndex == 2)
            % Short & Long trajectories
            t0       	=   ([0.31 0.32]-sfaParam.phi)/sfaParam.omega;
            
        else
            % Illegal Trajectory index -> Internal error
            error('sfaLP:sfaLPsolver:TrajType','Illegal trajectory type -- internal error')
        end
        
        % Find cut-off trajectory
        obj.t0_c        =   fminbnd(@(t0) -(sfaParam.potVector(tr(t0))-sfaParam.potVector(t0))^2,t0(1),t0(2));
        obj.tr_c        =   tr(obj.t0_c);
        obj.kappa       =   .5*(sfaParam.potVector(obj.tr_c)-sfaParam.potVector(obj.t0_c)).^2/sfaParam.Up;
        
        % Test for accuracy
        if abs(traj(obj.t0_c,obj.tr_c)) > 1e-10
            % Cut-off trajectory poorly determined
            warning('sfaLP:sfaLP:CutOffTraj','Cut-off trajectory poorly determined.')
        end
end
end
    end
    
%% Miscellaneous methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=protected)
function [t0,tr] = complexInitCond(obj,sfaParam,freq)
    % Approximate solution for complex trajectories beyond the cut-off
    %   with zero ionization potential (common to long and short traj.)
    
    % Scaled recollision energy
    DEr                 =   (freq-sfaParam.Ip)/sfaParam.Up-obj.kappa;
    
    % Approximate phases
    phi_0               =  -.0354*DEr -1i*.5400*acosh(.8057*DEr+1);
    phi_r               =   .0148*DEr +1i*.1895*acosh(.6916*DEr+1);
    
    % Convert to times
    t0                  =   (phi_0-sfaParam.phi)/sfaParam.omega+obj.t0_c;
    tr                  =   (phi_r-sfaParam.phi)/sfaParam.omega+obj.tr_c;
end 
    end
end

