classdef sfaLP_long < sfaLPsolver
%Linear polarization strong-field approximation solver for short
%trajectories. This class is part of the sfaLP suite and is not meant to be
%used on its own by end users.
    
    properties (Access=private)
    end
    
%% Inherited properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public,SetAccess=private)
        type                                % Trajectory type
    end
    
    properties (GetAccess=public,SetAccess=private,Hidden=true)
        t0_ref                              % Reference ionization time
        tr_ref                              % Reference recollision time
    end
    
    properties (Access=private)
        NLopts                              % Nonlinear solver options
    end
    
    
%% Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function obj = sfaLP_long(sfaParam)
% Constructor
    
    % Cut-off trajectory
    obj                 =   obj@sfaLPsolver(sfaParam,1);
    
    obj.type            =   'Long';
    obj.t0_m            =   -sfaParam.phi/sfaParam.omega;
    obj.tr_m            =   (2*pi-sfaParam.phi)/sfaParam.omega;
    
    % Reference values for future sfa computations
        % Nonlinear solver options
        obj.NLopts      =   optimset('Display','off','TolFun',1e-12);
        
        % Set reference frequencies
        freq          	=   sort(sfaParam.refFreq);

        if isempty(freq)
            freq        =   [.8*sfaParam.Ip 1.2*(obj.kappa*sfaParam.Up+sfaParam.Ip)];
        end

        if length(freq) == 2
            freq        =   linspace(freq(1),freq(2),sfaParam.refNbFreq);
        end
        
        % Compute reference frequency solutions
        if any(isnan(freq))
            % Do not set reference frequency for the problem
            obj.t0_ref  =   [];
            obj.tr_ref  =   [];
            
        else
            % Initialization
            [~,ind]     =   sort(abs(freq-2*sfaParam.Up-sfaParam.Ip));
            t0          =   NaN(size(freq));
            tr          =   NaN(size(freq));
            
            % Find first solution from no Ip
            if (freq(ind(1)) < sfaParam.Ip) || (freq(ind(1)) > obj.kappa*sfaParam.Up+sfaParam.Ip)
                error('sfaLP:sfaLP_long:RefFreq','Out-of-bound reference frequencies - does not overlap with the threshold to cut-off energy range');
            end
            
            tr(ind(1))  =   fzero(@(t) sfaParam.saddleEqRec(t,-sfaParam.potVector(obj.ionizationTimeReal(sfaParam,t)),freq(ind(1))),[obj.tr_m obj.tr_c]);
            t0(ind(1))	=   obj.ionizationTimeReal(sfaParam,tr(ind(1)));
            
            % Find first solution with accurate Ip
            X0          =   [t0(ind(1)) .001 tr(ind(1)) 0];
            
            for IpShift = .1:.1:sfaParam.Ip
                X0      =   obj.NLsolver(@(X) obj.complexSolution(sfaParam,X,freq(ind(1)),sfaParam.Ip-IpShift),X0,obj.NLopts);
            end
            
            X0          =   obj.NLsolver(@(X) obj.complexSolution(sfaParam,X,freq(ind(1)),0),X0,obj.NLopts);
            t0(ind(1))  =   X0(1) + 1i*X0(2);
            tr(ind(1))  =   X0(3) + 1i*X0(4);
            
            % Find following solutions
            for k=2:length(freq)
                % Get initial guess
                if k==2
                    % Use the previous guess
                    X0  =   [real(t0(ind(k-1))) imag(t0(ind(k-1))) real(tr(ind(k-1))) imag(tr(ind(k-1)))];
                elseif k==3
                    % Linear extrapolation fit
                    t_0 =   interp1(freq(ind(1:k-1)),t0(ind(1:k-1)),freq(k),'linear','extrap');
                    t_r =   interp1(freq(ind(1:k-1)),tr(ind(1:k-1)),freq(k),'linear','extrap');
                    X0  =   [real(t_0) imag(t_0) real(t_r) imag(t_r)];
                else
                    % Spline extrapolation fit
                    t_0 =   interp1(freq(ind(1:k-1)),t0(ind(1:k-1)),freq(k),'spline','extrap');
                    t_r =   interp1(freq(ind(1:k-1)),tr(ind(1:k-1)),freq(k),'spline','extrap');
                    X0  =   [real(t_0) imag(t_0) real(t_r) imag(t_r)];
                end
                
                % Compute solution
                X0      =   obj.NLsolver(@(X) obj.complexSolution(sfaParam,X,freq(ind(k)),0),X0,obj.NLopts);
                
                % Update data
                t0(ind(k))= X0(1) + 1i*X0(2);
                tr(ind(k))= X0(3) + 1i*X0(4);
                
            end
        end
        
        % Update properties
        if (length(freq) == 1) && ~isnan(freq)
            % Single value initial solution
            obj.t0_ref  =   @(f) t0;
            obj.tr_ref  =   @(f) tr;
            
        elseif length(freq) == 2
           % Linear fit
           obj.t0_ref   =   griddedInterpolant(freq,t0,'linear','linear');
           obj.tr_ref   =   griddedInterpolant(freq,tr,'linear','linear');
        elseif length(freq) > 2
            % Spline fit
           obj.t0_ref   =   griddedInterpolant(freq,t0,'spline','spline');
           obj.tr_ref   =   griddedInterpolant(freq,tr,'spline','spline');
        end
        
end
    end

%% Saddle point equation solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function [t0,tr,p] = saddlePointSolution (obj,sfaParam,Freq)
    % Solution of the saddle point equation
    
    % Initialization
    t0                  =   NaN(size(Freq));
    tr                  =   NaN(size(Freq));
    
    if isempty(obj.t0_ref)
        % Solver has not been properly initialized -> return error
        error('sfaLP:sfaLP_long:NoRefInit','Strong-field approximation solver has not been initialized with the appropriate reference(s)')
    end
    
    % Find solutions
    for k=1:numel(Freq)
        % Approximate solution from the reference
        t0(k)           =   obj.t0_ref(Freq(k));
        tr(k)           =   obj.tr_ref(Freq(k));
        X0              =   [real(t0(k)) imag(t0(k)) real(tr(k)) imag(tr(k))];
        
        % Find actual solution
        X0              =   obj.NLsolver(@(X) obj.complexSolution(sfaParam,X,Freq(k),0),X0,obj.NLopts);
                
        % Update data
        t0(k)           = X0(1) + 1i*X0(2);
        tr(k)           = X0(3) + 1i*X0(4);
    end
    
    % Compute momentum
    p                   =   sfaParam.saddlePointMomentum(t0,tr);
    
end
    end
    methods (Access=private)
function t0 = ionizationTimeReal(obj,sfaParam,tr)
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
function D = complexSolution(obj,sfaParam,X,freq,IpShift)                   %#ok<INUSL>
    % Compute sfa solution for complex trajectories
    
    % Initialization
    t0                  =   X(1) + 1i*X(2);
    tr                  =   X(3) + 1i*X(4);
    p                   =  	sfaParam.saddlePointMomentum(t0,tr);
    
    % Compute error
    dIon                =   sfaParam.saddleEqIon(t0,p) - IpShift;
    dRec                =   sfaParam.saddleEqRec(tr,p,freq);
    
    % Return result
    D                   =   [real(dIon);imag(dIon);real(dRec);imag(dRec)];
    
end
        
    end
    
end
