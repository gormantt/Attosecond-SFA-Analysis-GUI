classdef sfaLP < handle
%sfaLP Strong Field Approximation with Linear Polarization class 
%
%   This class models the strong-field approximation for linear
%   polarization, solving for the stationary phase equations
%       \nabla S(t_0,t_r,p) = 0                                         (1)
%   where the action is defined as
%       S(t_0,t_r,p) = -int( (p+A(s))^2/2 + Ip, s=t_0,t_r ) - nu t_r,   (2)
%   with int the integral, p the canonical momentum, A the vector 
%   potential, Ip the target ionization potential and nu the harmonic
%   frequency. In this formulation, t_0 and t_r are interpreted as the
%   ionization and recollision times, respectively.
%
%   From the action definition (2), the stationary phase equation (1) leads
%   to the system of equations
%       (p+A(t_0))^2/2 + Ip               = 0,                          (3)
%      -(p+A(t_r))^2/2 - Ip - nu          = 0,                          (4)
%      -p*(t_r-t_0) - int(A(s),s=t_0,t_r) = 0.                          (5)
%   Note that in order to satisfy Eq. (4) nu should be negative (the other
%   two terms begin positive). To simplify computations, we define the
%   input frequency as the opposite of it
%       freq = -nu.
%
%   In all equations, the linearly polarized laser is described by its
%   amplitude E0, frequency omega and carrier envelope phase phi, such that
%   the electric field/potential vector respectively read
%       E(t) =  E0 *       cos( omega*t + phi ),                        (6)
%       A(t) = -E0/omega * sin( omega*t + phi ).                        (7)
%   All parameters and variables in the class are expressed in atomic units
%   unless otherwise specified. Conversion to other sets of units can be
%   performed with with the cvUnits class suite of functions.
%
%sfaLP properties
%   LASER
%     * E0 - amplitude [ positive scalar {0.1} ]
%     * omega - frequency [ positive scalar {0.1} ]
%     * phi - carrier envelope phase [scalar {0} ]
%     # Up - ponderomotive energy [ positive scalar {N/A}]
%
%   TARGET
%     * Ip - ionization potential [ positive scalar {1} ]
%     # gamma - Keldysh parameter [ positive scalar {N/A} ]
%
%   SEMI-CLASSICAL TRAJECTORY
%     # kappa - maximum classical return energy in units of Up
%           [ positive scalar {N/A} ]
%     # trajType - trajectory type [ {Short} | Long | ShortNoIp | LongNoIp]
%       This parameter specifies the type of semi-classical trajectory to
%       search for in the stationary phase Eqs. (3-5).
%       'Short'         solves for the so-called short trajectories.
%       'Long'          solves for the so-called long trajectories.
%       'ShortNoIp'     solves for short trajectories while ignoring Ip
%                       (Ip=0) in Eq. (3), but still keeping it in (4).
%                       This leads to real-valued solution for frequencies
%                       between Ip and kappa*Up.
%       'ShortNoIp'     solves for long trajectories while ignoring Ip
%                       (Ip=0) in Eq. (3), but still keeping it in (4).
%                       This leads to real-valued solution for frequencies
%                       between Ip and kappa*Up.
%       Note that 'short' and 'long' trajectory, due to the complex nature
%       of the saddle point Eqs. (3-5), can have a big overhead at
%       initialization (when the class object is first created and
%       reinitialized). Depending on the use of its use, those can be
%       controled with the refFreq and refNbFreq properties (see below).
%
%   INITIALIZATION
%     * refFreq - reference frequencies [ {[]}, scalar, vector, NaN ]
%       This parameter is relevant only for 'short' and 'long'
%       trajectories. 
%       By nature, the stationary phase Eqs. (3-5) corresponds to a
%       nonlinear system of equation which we solve for numerically. Such
%       computations require for an initial guess. refFreq defines the
%       frequencies for which such initial guess should be computed and
%       stored for reference in future computations. Supported calues
%       []      let the program automatically define reference frequencies
%       scalar  computes and stores only one reference value - the one
%               defined by the scalar value - and uses this one for all
%               further computations.
%       [f1 f2] specifies the frequencies boundaries for the reference
%               computation. Reference frequencies are then defined as
%               linspace(f1,f2,refNbFreq) (see below).
%       vector  specifies the values of the frequencies to use for
%               reference.
%       NaN     (or if any component of the parameter is NaN) disables the
%               reference computation. This suppresses the overhead
%               complexity of the reference initialization. However, in
%               this case the methods saddlePointSolution and groupDelay
%               are not defined (see below) and will cause and error.
%       For practical and convergense consderations, with both the scalar
%       and vector options, at leas one of the reference frequencies mut be
%       between the classical threshold (Ip) and cut-off (kappa*Up+Ip).
%     * refNbFreq - number of reference frequencies
%           [ positive integer | {51} ]
%       For default and two-vector reference frequencies (refFreq = [] or
%       refFreq = [f1 f2] - see above), it defines the number of reference
%       frequencies to be used in the initialization.
%
%   GROUP DELAY
%     @ groupDelayStep - frequency step for group delay computations
%           [ positive scalar | {1e-4} ]
%       This parameter specifies the frequency step to be used for
%       (discrete derivative - centered, second order scheme) computation
%       of the group delay.
%
%   Properties marked with a '*' bullet have side effects on the definition
%   of the object and require a reinitialization after they are changed.
%   Their values can be accessed with the regular '.' operator. Changing
%   their values, however, requires the use of the 'set' method (see
%   below).
%
%   Properties marked with a '#' bullet are dependent and their values can
%   only be accessed and not changed (even with the 'set' method - see
%   below).
%
%   Properties marked with a '@' bullet can have their values both accessed
%   and changed with the regular '.' operator. They can also be changed
%   with the 'set' method but this will trigger a reinitialization of the
%   object (see below).
%
%sfaLP methods
%   INITIALIZATION
%       sfaLP('name1',value1,'name2',value2,...) creates a strong-field
%       approximation with linear polarization object in which the named
%       properties have the specified values. Any unspecified properties
%       have default values. Case is ignored for property. Valid names are
%       all the class properties with a '*' or '@' bullet plus 'type' to
%       define the trajectory type (same options as for the dependent
%       property trajType). Dependent properties ('#' bullet above) are not
%       valid names and cannot be set directly. This is the class
%       constructor.
%
%       sfaLP(optStruct) with an option structure optStruct as input uses
%       the values of the named fields as parameters. Valid names are the
%       same as above and are case sensitive. Additional fields to the
%       option structure are ignored.
%
%       sfaLP() or sfaLP([]) creates a strong-field approximation with
%       linear polarization object with all properties set to default
%       values.
%
%       set('name1',value1,'name2',value2,...) updates the named properties
%       with the specfied values. Any unspecified properties are left
%       unchanged. After updating the values, the object is reinitialized
%       to reflect changes (the cost of reinitialization is the same as for
%       object initialization with the given trajectory type). Valid names
%       follow the same rules as for the constructor with the corresponding
%       input variable pattern.
%
%       set(optStruct) with an option structure optStruct as input uses the
%       values of the named fields as parameters. This will also trigger a
%       reinitialization of the object. Field names to optStruct follow the
%       same rules as for the constructor with the corresponding input
%       variable pattern.
%
%   LASER FIELD
%       E = elecField(t) returns the value of the laser electric field at
%       time t, following Eq. (6).
%
%       A = potVector(t) returns the value of the laser potential vector at
%       time t, following Eq. (7).
%
%   SADDLE POINT EQUATIONS
%       SP = saddleEqIon(t0,p) returns the result of the ionization saddle
%       point Eq. (3) for ionization time t0 and momentum p. For solutions
%       of the saddle point system (3-5), the result should be 0 for
%       'short' and 'long' trajectories and Ip for their counterpart when
%       the ionization potential is ignored (real trajectories).
%
%       SP = saddleEqRec(tr,p,freq) returns the result of the recollision
%       saddle point Eq. (4) for recollision time tr, momentum p and
%       frequency freq. For solutions of the saddle point system (3-5), the
%       result should be 0.
%
%       SP = saddleEqTraj(t0,tr,p) returns the result of the trajectory
%       saddle point Eq. (5) for ionization and recollision times t0 and
%       tr, and momentum p. For solutions of the saddle point system (3-5),
%       the result should be 0 
%   
%   NOTES:
%     * sfaLP is a handle class (similar to plots). A handle is a reference
%       to an object. If you copy an object's handle variable, MATLAB 
%       copies only the handle. Both the original and copy refer to the
%       same object. For example, if a function modifies a handle object
%       passed as an input argument, the modification affects the original
%       input object.
%     * For 'short' and 'long' trajectory types, selecting the a scalar
%       reference frequency can significantly reduce the overhead
%       initialization cost. As a trade-off saddle point solution and group
%       delay computation might converge (to the accurate solution) only
%       around the reference.
%     * For 'short' and 'long' trajectories, selecting the vector reference
%       frequencies, it is advised to use a set that spans over all
%       frequencies that will be used in actual saddle point solution and
%       group delay computations. This is however not a requirement and the
%       program will extrapolate if required. During the initialization,
%       the program computes the reference solution by continuity across
%       the elements of refFreq. If the vector contains too wide gaps the
%       initialization might fail to converge (to the appropriate
%       solution), and in turn saddle point solutions and group delay
%       computations.
%     * Member methods saddleEqIon, saddleEqRec and saddleEqTraj can be
%       used to check the accuracy of computed solutions to the saddle
%       point system (3-5). 
%
%   See also cvUnits

% F. Mauger
%   Version 1.0.00
%   01/10/2017  Creation
%   01/11/2017  Clean up function and complete the documentation

%% System variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (GetAccess=public,SetAccess=private)
        % Laser
        %   E(t)	=	E0      *cos(omega*t+phi)
        %   A(t)    =  -E0/omega*sin(omega*t+phi)
        E0                                  % Peak laser field amplitude (a.u.)
        omega                               % Laser field frequency (a.u.)
        phi                                 % Laser field CEP  (a.u.)
        
        % Target
        Ip                                  % Ionization potential (a.u.)
        
        % Computations
        refFreq                             % Reference frequencies (a.u.)
        refNbFreq                           % Default number of reference frequencies
    end
    
    properties (Access=public)
        % Group delay
        groupDelayStep                      % Energy step for group delay computation (a.u.)
    end
    
    properties (Dependent)
        Up                                  % Ponderomotive energy (dependent - a.u.)
        gamma                               % Keldysh parameter (dependent - a.u.)
        trajType                            % Trajectory type
        kappa                               % Maximum return energy (dependent - Up)
    end
    
    properties (Access=private)
        solver                              % The actual solver
        type                                % A copy of the type name, to handle type update
    end
    
%% Dependent properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        function Up = get.Up(obj),      Up = (.5*obj.E0/obj.omega)^2; end
        function ga = get.gamma(obj),   ga = sqrt(.5*obj.Ip/obj.Up); end
        function tT = get.trajType(obj),tT = obj.solver.type; end
        function ka = get.kappa(obj),   ka = obj.solver.kappa; end
    end
    
%% Class functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function obj = sfaLP(varargin)
% sfaLP strong-field approximation with linear polarization constructor
%
%       sfaLP('name1',value1,'name2',value2,...) creates a strong-field
%       approximation with linear polarization object in which the named
%       properties have the specified values. Any unspecified properties
%       have default values. Case is ignored for property. Valid names
%       correspond to all class member propertie (except dependent ones)
%       plus 'type' to specify the trajectory type. See the class
%       documentation for a list of valid names and their description. 
%
%       sfaLP(optStruct) with an option structure optStruct as input uses
%       the values of the named fields as parameters. Valid names are the
%       same as above and are case sensitive. Additional fields to the
%       option structure are ignored.
    
    % Set default variables
    obj.E0              =   .1;
    obj.omega           =   .1;
    obj.phi             =   0;
    
    obj.Ip              =   1;
    
    obj.refFreq         =   [];
    obj.refNbFreq       =   51;
    
    obj.groupDelayStep  =   1e-4;
    
    obj.type            =   'Short';
    
    % Set parameters to their actual values
    %   missing parameters will be generated in there
    obj.set(varargin{:});
    
end
function set(obj,varargin)
% set set givenmember properties to selected values
%
%       set('name1',value1,'name2',value2,...) updates the named properties
%       with the specfied values. Any unspecified properties are left
%       unchanged. After updating the values, the object is reinitialized
%       to reflect changes (the cost of reinitialization is the same as for
%       object initialization with the given trajectory type). Valid names
%       correspond to all class member propertie (except dependent ones)
%       plus 'type' to specify the trajectory type. See the class
%       documentation for a list of valid names and their description. 
%       input variable pattern.
%
%       set(optStruct) with an option structure optStruct as input uses the
%       values of the named fields as parameters. This will also trigger a
%       reinitialization of the object. Valid names are the same as above
%       and are case sensitive. Additional fields to the option structure
%       are ignored.

%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization
    Names               =   properties(obj);
    Names               =   setdiff(Names,{'Up','gamma','trajType','kappa'});
    Names               =   [Names;{'type'}];
    
    % Update relevant fields
    if (nargin==2) && (isstruct(varargin{1}))
        % Structure input
        in              =   varargin{1};
        
        for k=1:length(Names)
            if isfield(in,Names{k})
                obj.(Names{k}) = in.(Names{k});
            end
        end
        
    elseif mod(nargin,2) == 1
        % Name-value pairs input
        for k=2:2:nargin
            % Test for option type
            if ~ischar(varargin{k-1})
                error('sfaLP:sfaLPset:NoPropName',['Property name expected in argument # ' num2str(k-1) '.']);
            end

            % Update appropriate field
            is_field	=   false;
            for l=1:length(Names)
                if strcmpi(varargin{k-1},Names{l})
                    % We have a match -> update data and proceed to next pair
                    obj.(Names{l}) = varargin{k};
                    is_field=true;
                    break;
                end
            end

            if ~is_field
                warning('sfaLP:sfaLP:InvalidPropName',['Invalid property name ' varargin{k-1} '. Property ignored']);
            end
        end
        
    else
        % Return error
        error('sfaLP:sfaLP:InputArg','Input arguments must be a structure of parameters or occur in name-value pairs.');
    end
    
    % Update solver (and type)
    switch lower(obj.type)
        case lower({'Short, no Ip','ShortNoIp','short_no_Ip'})
            % Short trajectories, ignoring Ip
            obj.solver  =   sfaLP_shortNoIp(obj);
            
        case lower({'Long, no Ip','LongNoIp','long_no_Ip'})
            % Short trajectories, ignoring Ip
            obj.solver  =   sfaLP_longNoIp(obj);
            
        case lower({'short','S'})
            % Short trajectories
            obj.solver  =   sfaLP_short(obj);
            
        case lower({'long','L'})
            % Short trajectories
            obj.solver  =   sfaLP_long(obj);
            
        otherwise
            error('sfaLP:sfaLP:InvalidTrajType',['Invalid trajectory type ' opts.type '.']);
    end
    
    obj.type            =   obj.solver.type;
    
end
    end

%% Laser field functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function E = elecField(obj,t)
%elecField laser electric field
%
%   E = elecField(t) returns the value of the laser electric field at given 
%   time t
%       E(t) = E0*cos(omega*t+phi).
    E                   =   obj.E0*cos(obj.omega*t+obj.phi);
end
function A = potVector(obj,t)
%potVector laser potential vector
%
%   A = potVector(t) returns the value of the laser potential vector at
%   given time t
%       A(t) = -E0/omega*sin(omega*t+phi).
    A                   =   -obj.E0/obj.omega*sin(obj.omega*t+obj.phi);
end
    end

%% Saddle Point equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function  SP = saddleEqIon(obj,t0,p)
    % Saddle point equation for ionization (3). For trajectories solutions
    % of the saddle point approximation, the result should be 0 or Ip if
    % the ionization potential is ignored in the solution.
    SP                  =   .5*(p+obj.potVector(t0)).^2 + obj.Ip;
end
function SP = saddleEqRec(obj,tr,p,freq)
    % Saddle point equation for recollision (4). For trajectories solutions
    % of the saddle point approximation, the result should be 0.
    SP                  =   .5*(p+obj.potVector(tr)).^2 + obj.Ip - freq;
end
function SP = saddleEqTraj(obj,t0,tr,p)
    % Saddle point equation for trajectory (5). For trajectories solutions
    % of the saddle point approximation, the result should be 0.
    SP                  =   p.*(tr-t0) + obj.E0/obj.omega^2*( cos(obj.omega*tr+obj.phi) - cos(obj.omega*t0+obj.phi));
end
    end
    
%% Saddle point solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function [t0,tr,p,S] = saddlePointSolution(obj,freq)
    % Compute solutions of the saddle point equations for given harmonic
    % frequencies.
    
    % Find solutions with appropriate solver
    [t0,tr,p]           =   obj.solver.saddlePointSolution(obj,freq);
    
    % Compute SFA phases
    [Sa,Sc,Sf]          =   obj.saddlePointPhase(t0,tr,p,freq);
    S                   =   Sa+Sc+Sf;
end
function [Sact,Score,Sfreq] = saddlePointPhase(obj,t0,tr,p,Freq)
    % Compute the phase components of the strong-field approximation saddle
    % point solution, respectively defined as
    %   action          -int((p+A(s))^2/2,s=t0,tr),
    %   core            -Ip*(tr-t0),
    %   frequency       -nu*tr,
    % with nu = -freq. The total phase is then defined as the sum of the
    % three aforementioned contributions.
    %
    % By definition, theintegrand of the  action component is
    %   .5*(p-E0/omega*sin(omega*t+phi))^2,
    %       = .5*p^2 - p*E0/omega*sin(omega*t+phi)
    %         	+ .5*E0^2/omega^2*sin(omega*t+phi)^2,
    %       = .5*p^2 + .25*E0^2/omega^2 - p*E0/omega*sin(omega*t+phi)
    %           -.25*E0^2/omega^2*cos(2*omega*t+2*phi).
    % Now integrating between ionization and recollision times, we get the
    % action component
    %   Sact= -(.5*p^2+.25*E0^2/omega^2)*(tr-t0)
    %           - p*E0/omega^2*(cos(omega*tr+phi)-cos(omega*t0+phi))
    %           +.125*E0^2/omega^3*(sin(2*omega*tr+2*phi)-sin(2*omega*t0+2*phi))
    Sact                =   -(.5*p.^2+.25*obj.E0^2/obj.omega^2).*(tr-t0) ...
                            -obj.E0/obj.omega^2*p.*(cos(obj.omega*tr+obj.phi)-cos(obj.omega*t0+obj.phi)) ...
                            +.125*obj.E0^2/obj.omega^3*(sin(2*obj.omega*tr+2*obj.phi)-sin(2*obj.omega*t0+2*obj.phi));
    Score               =  -obj.Ip*(tr-t0);
    Sfreq               =   Freq.*tr;
end
function p = saddlePointMomentum(obj,t0,tr)
    % Compute the momentum to satisfy the saddle point recollision Eq. (5),
    % such that 
    %   p   =   -int(A(s),s=t_0,t_r)/(t_r-t_0),
    %           -E0/omega^2*(cos(omega*t_r+phi)-cos(omega*t_0+phi))/(t_r-t_0).
    p                   =   -obj.E0/obj.omega^2*( cos(obj.omega*tr+obj.phi) - cos(obj.omega*t0+obj.phi))./(tr-t0);
end
    end
%% Group delay computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
function GD = groupDelay(obj,Freq)
    % Compute the group delay at given frequency using a second order
    % discrete derivative scheme.
    
    % Group delay computation
    [~,~,~,S_down]      =   obj.saddlePointSolution(Freq-obj.groupDelayStep);
    [~,~,~,S_up]        =   obj.saddlePointSolution(Freq+obj.groupDelayStep);
    
    GD                  =   .5*real(S_up-S_down)/obj.groupDelayStep;

end
    end
    
end
