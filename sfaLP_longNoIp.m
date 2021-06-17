classdef sfaLP_longNoIp < sfaLPsolver
%Linear polarization strong-field approximation solver for short
%trajectories, ignoring the ionization potential in the ionization step. 
%This class is part of the sfaLP suite and is not meant to be used on its
%own by end users.
    
    properties (Access=private)
    end
    
%% Inherited properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public,SetAccess=private)
        type                                % Trajectory type
    end
    
%% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function obj = sfaLP_longNoIp(sfaParam)
% Constructor
    
    % Cut-off trajectory
    obj                 =   obj@sfaLPsolver(sfaParam,1);
    
    obj.type            =   'LongNoIp';
    obj.t0_m            =   -sfaParam.phi/sfaParam.omega;
    obj.tr_m            =   (2*pi-sfaParam.phi)/sfaParam.omega;
    
    % Nothing else to do here
end
    end
    
%% Saddle point equation solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function [t0,tr,p] = saddlePointSolution (obj,sfaParam,Freq)
    % Solution of the saddle point equation
    
    % Initialization
    t0                  =   NaN(size(Freq));
    tr                  =   NaN(size(Freq));
    
    % Compute sfa trajectories
    for k=1:numel(Freq)
        if Freq(k) < sfaParam.Ip
            % Below threshold trajectories are not defined
            warning('sfaLP:sfaLP_longNoIp:BelowThreshold','Below threshold trajectories are not defined, NaN returned instead');
        elseif Freq(k) < sfaParam.Ip + sfaParam.Up*1e-11
            % Zero return energy -> we know the solution
            t0(k)       =   obj.t0_m;
            tr(k)       =   obj.tr_m;
        elseif Freq(k) < (1-1e-11)*obj.kappa*sfaParam.Up+sfaParam.Ip
            % Real valued trajectory solutions
            tr(k)       =   fzero(@(t) sfaParam.saddleEqRec(t,-sfaParam.potVector(obj.ionizationTime(sfaParam,t)),Freq(k)),[obj.tr_m obj.tr_c]);
            t0(k)       =   obj.ionizationTime(sfaParam,tr(k));
        elseif Freq(k) < (1+1e-11)*obj.kappa*sfaParam.Up+sfaParam.Ip
            % Cut-off trajectory -> we know the solution
            t0(k)       =   obj.t0_c;
            tr(k)       =   obj.tr_c;
        else
            % Above cut-off trajectories -> complex solutions
            [t0(k),tr(k)]=  obj.complexInitCond(sfaParam,Freq(k));
            
            X0          =   [real(t0(k));imag(t0(k));real(tr(k));imag(tr(k))];
            X           =   obj.NLsolver(@(X) obj.complexSolution(sfaParam,X,Freq(k)),X0,optimset('Display','off','TolFun',1e-10));
            
            t0(k)       =   X(1) + 1i*X(2);
            tr(k)       =   X(3) + 1i*X(4);
        end
        
    end
    
    % Compute momentum
    p                   =  -sfaParam.potVector(t0);
    
end
    end
    methods (Access=private)
function t0 = ionizationTime(obj,sfaParam,tr)
    % Computes ionization time for a given recollision time
    
    if tr < obj.tr_c + 1e-11
        % Threshold trajectory -> We know the ionization time
        t0              =   obj.t0_c;
    elseif tr > obj.tr_m - 1e-11
        % Cut-off trajectory -> We know the recollision time
        t0              =   obj.t0_m;
    else
        % Find trajectory
        t0              =   fzero(@(t0) sfaParam.saddleEqTraj(t0,tr,-sfaParam.potVector(t0)),[obj.t0_m obj.t0_c]);
    end
end
function D = complexSolution(obj,sfaParam,X,freq)                           %#ok<INUSL>
        % Compute sfa solution for complex trajectories
        
        % Initialization
        t0              =   X(1) + 1i*X(2);
        tr              =   X(3) + 1i*X(4);
        p               =  -sfaParam.potVector(t0);
        
        % Compute errors
        dRec            =   sfaParam.saddleEqRec(tr,p,freq);
        dTraj           =   sfaParam.saddleEqTraj(t0,tr,p);
        
        % Return result
        D               =   [real(dRec);imag(dRec);real(dTraj);imag(dTraj)];
end
    end
    
end

