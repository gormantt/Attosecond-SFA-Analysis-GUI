classdef systemSolve
% Wrapper for nonliner solver (Matlab fsolve)
    
    methods (Static=true)
function [x,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUN,x,options,varargin)
%FSOLVE solves systems of nonlinear equations of several variables.
%
%   FSOLVE attempts to solve equations of the form:
%             
%   F(X)=0    where F and X may be vectors or matrices.   
%
%   X=FSOLVE(FUN,X0) starts at the matrix X0 and tries to solve the 
%   equations in FUN.  FUN accepts input X and returns a vector (matrix) of 
%   equation values F evaluated at X. 
%
%   X=FSOLVE(FUN,X0,OPTIONS) solves the equations with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%   FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType,
%   NonlEqnAlgorithm, MaxFunEvals, MaxIter, PlotFcns, OutputFcn,
%   DiffMinChange and DiffMaxChange, LargeScale, MaxPCGIter,
%   PrecondBandWidth, TolPCG, and TypicalX. Use the Jacobian option to
%   specify that FUN also returns a second output argument J that is the
%   Jacobian matrix at the point X. If FUN returns a vector F of m
%   components when X has length n, then J is an m-by-n matrix where J(i,j)
%   is the partial derivative of F(i) with respect to x(j). (Note that the
%   Jacobian J is the transpose of the gradient of F.)
%
%   X = FSOLVE(PROBLEM) solves system defined in PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fsolve' in PROBLEM.solver.  Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. The structure PROBLEM 
%   must have all the fields.
%
%   [X,FVAL]=FSOLVE(FUN,X0,...) returns the value of the equations FUN at X. 
%
%   [X,FVAL,EXITFLAG]=FSOLVE(FUN,X0,...) returns an EXITFLAG that describes the
%   exit condition of FSOLVE. Possible values of EXITFLAG and the corresponding 
%   exit conditions are
%
%     1  FSOLVE converged to a solution X.
%     2  Change in X smaller than the specified tolerance.
%     3  Change in the residual smaller than the specified tolerance.
%     4  Magnitude of search direction smaller than the specified tolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Algorithm terminated by the output function.
%    -2  Algorithm seems to be converging to a point that is not a root.
%    -3  Trust region radius became too small.
%    -4  Line search cannot sufficiently decrease the residual along the current
%         search direction.
%
%   [X,FVAL,EXITFLAG,OUTPUT]=FSOLVE(FUN,X0,...) returns a structure OUTPUT
%   with the number of iterations taken in OUTPUT.iterations, the number of
%   function evaluations in OUTPUT.funcCount, the algorithm used in OUTPUT.algorithm,
%   the number of CG iterations (if used) in OUTPUT.cgiterations, the first-order 
%   optimality (if used) in OUTPUT.firstorderopt, and the exit message in
%   OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,JACOB]=FSOLVE(FUN,X0,...) returns the 
%   Jacobian of FUN at X.  
%
%   Examples
%     FUN can be specified using @:
%        x = fsolve(@myfun,[2 3 4],optimset('Display','iter'))
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x);
%
%   FUN can also be an anonymous function:
%
%       x = fsolve(@(x) sin(3*x),[1 4],optimset('Display','off'))
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the system of 
%   nonlinear equations given in the function myfun, which is parameterized 
%   by its second argument c. Here myfun is an M-file function such as
%     
%       function F = myfun(x,c)
%       F = [ 2*x(1) - x(2) - exp(c*x(1))
%             -x(1) + 2*x(2) - exp(c*x(2))];
%           
%   To solve the system of equations for a specific value of c, first assign the
%   value to c. Then create a one-argument anonymous function that captures 
%   that value of c and calls myfun with two arguments. Finally, pass this anonymous
%   function to FSOLVE:
%
%       c = -1; % define parameter first
%       x = fsolve(@(x) myfun(x,c),[-5;-5])
%
%   See also OPTIMSET, LSQNONLIN, @, INLINE.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:18 $

% ------------Initialization----------------

defaultopt = struct('Display','final','LargeScale','off',...
   'NonlEqnAlgorithm','dogleg',...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off','FunValCheck','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','OutputFcn',[],'PlotFcns',[]);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 3, options=[]; end

% Detect problem structure input
if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error('optim:fsolve:InputArg','The input to FSOLVE should be either a structure with valid fields or consist of at least two arguments.');
    end
end

if nargin == 0
  error('optim:fsolve:NotEnoughInputs','FSOLVE requires at least two input arguments.')
end

% Check for non-double inputs
if ~isa(x,'double')
  error('optim:fsolve:NonDoubleInput', ...
        'FSOLVE only accepts inputs of data type double.')
end

LB = []; UB = []; 
xstart=x(:);
numberOfVariables=length(xstart);

large        = 'large-scale';
medium       = 'medium-scale: line search';
dogleg       = 'trust-region dogleg';

switch optimget(options,'Display',defaultopt,'fast')
    case {'off','none'}
        verbosity = 0;
    case 'iter'
        verbosity = 2;
    case 'final'
        verbosity = 1;
    case 'testing'
        verbosity = Inf;
    otherwise
        verbosity = 1;
end
diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
% 0 means large-scale trust-region, 1 means medium-scale algorithm
mediumflag = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
switch optimget(options,'NonlEqnAlgorithm',defaultopt,'fast')
    case 'dogleg'
        algorithmflag = 1;
    case 'lm'
        algorithmflag = 2;
    case 'gn'
        algorithmflag = 3;
    otherwise
        algorithmflag = 1;
end
mtxmpy = optimget(options,'JacobMult',defaultopt,'fast');
if isequal(mtxmpy,'atamult')
    warning('optim:fsolve:NameClash', ...
        ['Potential function name clash with a Toolbox helper function:\n' ...
        'Use a name besides ''atamult'' for your JacobMult function to\n' ...
        'avoid errors or unexpected results.'])
end

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    funfcn = systemSolve.lsqfcnchk(FUN,'fsolve',length(varargin),funValCheck,gradflag);
else
    error('optim:fsolve:InvalidFUN', ...
        ['FUN must be a function name, valid string expression, or inline object;\n' ...
        ' or, FUN may be a cell array that contains these type of objects.'])
end

JAC = [];
x(:) = xstart;
switch funfcn{1}
    case 'fun'
        fuser = feval(funfcn{3},x,varargin{:});
        f = fuser(:);
        nfun=length(f);
    case 'fungrad'
        [fuser,JAC] = feval(funfcn{3},x,varargin{:});
        f = fuser(:);
        nfun=length(f);
    case 'fun_then_grad'
        fuser = feval(funfcn{3},x,varargin{:});
        f = fuser(:);
        JAC = feval(funfcn{4},x,varargin{:});
        nfun=length(f);
    otherwise
        error('optim:fsolve:UndefinedCalltype','Undefined calltype in FSOLVE.')
end

if gradflag
    % check size of JAC
    [Jrows, Jcols]=size(JAC);
    if isempty(mtxmpy)
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows~=nfun || Jcols~=numberOfVariables
            error('optim:fsolve:InvalidJacobian', ...
                ['User-defined Jacobian is not the correct size:\n' ...
                ' the Jacobian matrix should be %d-by-%d.'],nfun,numberOfVariables)
        end
    end
else
    Jrows = nfun;
    Jcols = numberOfVariables;
end

XDATA = []; YDATA = []; caller = 'fsolve';

% Choose what algorithm to run: determine (i) OUTPUT.algorithm and 
% (ii) if and only if OUTPUT.algorithm = medium, also option.LevenbergMarquardt.
% Option LevenbergMarquardt is used internally; it's not user settable. For
% this reason we change this option directly, for speed; users should use
% optimset.
if ~mediumflag 
    if nfun >= numberOfVariables
        % large-scale method and enough equations (as many as variables)
        OUTPUT.algorithm = large;
    else 
        % large-scale method and not enough equations - switch to medium-scale algorithm
        warning('optim:fsolve:FewerFunsThanVars', ...
                ['Large-scale method requires at least as many equations as variables;\n' ...
                 ' using line-search method instead.'])
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'off'; 
    end
else
    if algorithmflag == 1 && nfun == numberOfVariables 
        OUTPUT.algorithm = dogleg;
    elseif algorithmflag == 1 && nfun ~= numberOfVariables
        warning('optim:fsolve:NonSquareSystem', ...
                ['Default trust-region dogleg method of FSOLVE cannot\n handle non-square systems; ', ...
                 'using Gauss-Newton method instead.']);
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'off';
    elseif algorithmflag == 2
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'on';
    else % algorithmflag == 3
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'off';
    end
end

if diagnostics > 0
    % Do diagnostics on information so far
    constflag = 0; gradconstflag = 0; non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
    confcn{1}=[];c=[];ceq=[];cGRAD=[];ceqGRAD=[];
    hessflag = 0; HESS=[];
    diagnose('fsolve',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
        mediumflag,options,defaultopt,xstart,non_eq,...
        non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn,f,JAC,HESS,c,ceq,cGRAD,ceqGRAD);

end

% Execute algorithm
if isequal(OUTPUT.algorithm, large)
    if ~gradflag
        Jstr = optimget(options,'JacobPattern',defaultopt,'fast');
        if ischar(Jstr)
            if isequal(lower(Jstr),'sparse(ones(jrows,jcols))')
                Jstr = sparse(ones(Jrows,Jcols));
            else
                error('optim:fsolve:InvalidJacobPattern', ...
                    'Option ''JacobPattern'' must be a matrix if not the default.')
            end
        end
    else
        Jstr = [];
    end
    computeLambda = 0;
    [x,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msg]=...
        snls(funfcn,x,LB,UB,verbosity,options,defaultopt,f,JAC,XDATA,YDATA,caller,...
        Jstr,computeLambda,varargin{:});
elseif isequal(OUTPUT.algorithm, dogleg)
    % trust region dogleg method
    Jstr = [];
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg]=...
        systemSolve.trustnleqn(funfcn,x,verbosity,gradflag,options,defaultopt,f,JAC,...
        Jstr,varargin{:});
else
    % line search (Gauss-Newton or Levenberg-Marquardt)
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg] = ...
        systemSolve.nlsq(funfcn,x,verbosity,options,defaultopt,f,JAC,XDATA,YDATA,caller,varargin{:});
end

Resnorm = FVAL'*FVAL;  % assumes FVAL still a vector
if EXITFLAG > 0 % if we think we converged:
    if Resnorm > sqrt(optimget(options,'TolFun',defaultopt,'fast'))
        OUTPUT.message = ...
            sprintf(['Optimizer appears to be converging to a minimum that is not a root:\n' ...
            'Sum of squares of the function values is > sqrt(options.TolFun).\n' ...
            'Try again with a new starting point.']);
        if verbosity > 0
            disp(OUTPUT.message)
        end
        EXITFLAG = -2;
    else
        OUTPUT.message = msg;
        if verbosity > 0
            disp(OUTPUT.message);
        end
    end
else
    OUTPUT.message = msg;
    if verbosity > 0
        disp(OUTPUT.message);
    end
end

% Reset FVAL to shape of the user-function output, fuser
FVAL = reshape(FVAL,size(fuser));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [allfcns,msg] = lsqfcnchk(funstr,caller,lenVarIn,funValCheck,gradflag)
%LSQFCNCHK Pre- and post-process function expression for FCNCHK.
%   [ALLFCNS,MSG] = LSQFUNCHK(FUNSTR,CALLER,lenVarIn,GRADFLAG) takes
%   the (nonempty) expression FUNSTR from CALLER with LenVarIn extra arguments,
%   parses it according to what CALLER is, then returns a string or inline
%   object in ALLFCNS.  If an error occurs, this message is put in MSG.
%
%   ALLFCNS is a cell array:
%    ALLFCNS{1} contains a flag
%    that says if the objective and gradients are together in one function
%    (calltype=='fungrad') or in two functions (calltype='fun_then_grad')
%    or there is no gradient (calltype=='fun'), etc.
%    ALLFCNS{2} contains the string CALLER.
%    ALLFCNS{3}  contains the objective function
%    ALLFCNS{4}  contains the gradient function (transpose of Jacobian).
%
%    If funValCheck is 'on', then we update the funfcn's (fun/Jacobian) so
%    they are called through CHECKFUN to check for NaN's, Inf's, or complex
%    values. Add a wrapper function, CHECKFUN, to check for NaN/complex
%    values without having to change the calls that look like this:
%    f = funfcn(x,varargin{:}); x is the first argument to CHECKFUN, then
%    the user's function, then the elements of varargin. To accomplish this
%    we need to add the user's function to the beginning of varargin, and
%    change funfcn to be CHECKFUN.
%

%   Copyright 1990-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:47:39 $

%    NOTE: we assume FUNSTR is nonempty.
% Initialize
msg='';
allfcns = {};
funfcn = [];
gradfcn = [];

if gradflag
    calltype = 'fungrad';
else
    calltype = 'fun';
end

% {fun}
if isa(funstr, 'cell') && length(funstr)==1
    % take the cellarray apart: we know it is nonempty
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, idandmsg] = fcnchk(funstr{1},lenVarIn);
    % Insert call to nested function checkfun which calls user funfcn
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end

    if ~isempty(idandmsg)
        error(idandmsg)
    end

    % {fun,[]}
elseif isa(funstr, 'cell') && length(funstr)==2 && isempty(funstr{2})
    if gradflag
        calltype = 'fungrad';
    end
    [funfcn, idandmsg] = fcnchk(funstr{1},lenVarIn);
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end
    if ~isempty(idandmsg)
        error(idandmsg)
    end

    % {fun, grad}
elseif isa(funstr, 'cell') && length(funstr)==2 % and ~isempty(funstr{2})

    [funfcn, idandmsg] = fcnchk(funstr{1},lenVarIn);
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end
    if ~isempty(idandmsg)
        error(idandmsg)
    end
    [gradfcn, idandmsg] = fcnchk(funstr{2},lenVarIn);
    if funValCheck
        userfcn = gradfcn;
        gradfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end

    if ~isempty(idandmsg)
        error(idandmsg)
    end
    calltype = 'fun_then_grad';
    if ~gradflag
        warning('optim:lsqfcnchk:IgnoringJacobian', ...
            ['Jacobian function provided but OPTIONS.Jacobian=''off'';\n' ...
            ' ignoring Jacobian function and using finite-differencing.\n' ...
            ' Rerun with OPTIONS.Jacobian=''on'' to use Jacobian function.'])
        calltype = 'fun';
    end

elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string or inline object
    [funfcn, idandmsg] = fcnchk(funstr,lenVarIn);
    if funValCheck
        userfcn = funfcn;
        funfcn = @checkfun; %caller and userfcn are in scope in nested checkfun
    end

    if ~isempty(idandmsg)
        error(idandmsg)
    end
    if gradflag % gradient and function in one function/M-file
        gradfcn = funfcn; % Do this so graderr will print the correct name
    end
else
    error('optim:lsqfcnchk:InvalidFUN', ...
        ['FUN must be a function or an inline object;\n', ...
        ' or, FUN may be a cell array that contains these type of objects.']);
end

allfcns{1} = calltype;
allfcns{2} = caller;
allfcns{3} = funfcn;
allfcns{4} = gradfcn;
allfcns{5}=[];

%------------------------------------------------------------
    function [f,J] = checkfun(x,varargin)
        % CHECKFUN checks for complex or NaN results from userfcn.
        % Inputs CALLER and USERFCN come from scope of OPTIMFCNCHK.
        % We do not make assumptions about f, or J. For generality, assume
        % they can all be matrices.

        if nargout == 1
            f = userfcn(x,varargin{:});
            if any(any(isnan(f)))
                error('optim:lsqfcnchk:checkfun:NaNFval', ...
                    'User function ''%s'' returned NaN when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif ~isreal(f)
                error('optim:lsqfcnchk:checkfun:ComplexFval', ...
                    'User function ''%s'' returned a complex value when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif any(any(isinf(f)))
                error('optim:lsqfcnchk:checkfun:InfFval', ...
                    'User function ''%s'' returned Inf or -Inf when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            end

        elseif nargout == 2
            [f,J] = userfcn(x,varargin{:});
            if any(any(isnan(f))) || any(any(isnan(J)))
                error('optim:lsqfcnchk:checkfun:NaNFval', ...
                    'User function ''%s'' returned NaN when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif ~isreal(f) || ~isreal(J)
                error('optim:lsqfcnchk:checkfun:ComplexFval', ...
                    'User function ''%s'' returned a complex value when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            elseif any(any(isinf(f))) || any(any(isinf(J)))
                error('optim:lsqfcnchk:checkfun:InfFval', ...
                    'User function ''%s'' returned Inf or -Inf when evaluated;\n %s cannot continue.', ...
                    functiontostring(userfcn),upper(caller));
            end
        end

    end %checkfun

end % lsqfcnchk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,Fvec,JAC,EXITFLAG,OUTPUT,msg]=...
  trustnleqn(funfcn,x,verbosity,gradflag,options,defaultopt,Fvec,JAC,...
       JACfindiff,varargin)
%TRUSTNLEQN Trust-region dogleg nonlinear systems of equation solver.
%
%   TRUSTNLEQN solves a system of nonlinear equations using a dogleg trust
%   region approach.  The algorithm implemented is similar in nature
%   to the FORTRAN program HYBRD1 of J.J. More', B.S.Garbow and K.E. 
%   Hillstrom, User Guide for MINPACK 1, Argonne National Laboratory, 
%   Rept. ANL-80-74, 1980, which itself was based on the program CALFUN 
%   of M.J.D. Powell, A Fortran subroutine for solving systems of
%   nonlinear algebraic equations, Chap. 7 in P. Rabinowitz, ed.,
%   Numerical Methods for Nonlinear Algebraic Equations, Gordon and
%   Breach, New York, 1970.

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/12/15 19:29:13 $
%   Richard Waltz, June 2001
%
% NOTE: 'x' passed in and returned in matrix form.
%       'Fvec' passed in and returned in vector form.
%
% Throughout this routine 'x' and 'F' are matrices while
% 'xvec', 'xTrial', 'Fvec' and 'FTrial' are vectors. 
% This was done for compatibility with the 'fsolve.m' interface.

% Define some sizes.
xvec = x(:);         % vector representation of x
nfnc = length(Fvec);  
nvar = length(xvec);

% Get user-defined options.
[maxfunc,maxit,tolf,tolx,derivCheck,DiffMinChange,...
     DiffMaxChange,mtxmpy,typx,giventypx,JACfindiff,structure,outputfcn,plotfcns] = ...
     systemSolve.getOpts(nfnc,nvar,options,defaultopt,gradflag);
if giventypx    % scaling featured only enabled when typx values provided
  scale = true;    
else
  scale = false;
end
broyden = false;    % broyden feature will be chosen by user in the future

% Handle the output function
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end
stop = false;

% Handle the plot function
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

% Initialize local arrays.
d       = zeros(nvar,1);
scalMat = ones(nvar,1); 
grad    = zeros(nvar,1);
JACd    = zeros(nvar,1);
xTrial  = zeros(nvar,1);
F       = zeros(size(x));
FTrial  = zeros(nvar,1);
if derivCheck
    if gradflag, 
        JACfindiff = JAC; % Initialize finite difference Jacobian with 
    else                % structure given by real Jacobian 
        if verbosity > 0              
            warning('optim:trustnleqn:DerivativeCheckOff', ...
                    ['DerivativeCheck on but analytic Jacobian not provided;\n' ...
                     '         turning DerivativeCheck off.'])
        end
        derivCheck = false;
    end
end

% Initialize some trust region parameters.
Delta    = 1e0;
DeltaMax = 1e10;
eta1     = 0.05;
eta2     = 0.9;
alpha1   = 2.5;
alpha2   = 0.25;

% Other initializations.
iter = 0;
numFevals = 1;   % computed in fsolve.m
numFDfevals = 0;
if gradflag
  numJevals = 1; % computed in fsolve.m
else
  numJevals = 0;
end 
done = false;
stepAccept = true;
numReject = 0;
exitStatus = 0;
normd = 0.0e0;
normdscal = 0.0e0;
scalemin = eps;
scalemax = 1/scalemin;
objold = 1.0e0;
broydenJac = false;
broydenStep = false;
obj = 0.5*Fvec'*Fvec;  % Initial Fvec computed in fsolve.m

% Compute initial finite difference Jacobian, objective and gradient.
if derivCheck || ~gradflag
  if structure && issparse(JACfindiff)
    group = color(JACfindiff); % only do color if given some structure and sparse
  else
    group = 1:nvar;
  end 
  [JACfindiff,numFDfevals] = systemSolve.sfdnls(x,Fvec,JACfindiff,group,[], ...
                      DiffMinChange,DiffMaxChange,funfcn{3},[],[],varargin{:});
  numFevals = numFevals + numFDfevals;
end

switch funfcn{1}
case 'fun'
  JAC = JACfindiff;
case 'fungrad'         % Initial Jacobian computed in fsolve.m
  if derivCheck, graderr(JACfindiff,JAC,funfcn{3}); end
case 'fun_then_grad'   % Initial Jacobian computed in fsolve.m
  if derivCheck, graderr(JACfindiff,JAC,funfcn{4}); end
otherwise
  error('optim:trustnleqn:UndefinedCalltype','Undefined calltype in FSOLVE.')
end 
grad = feval(mtxmpy,JAC,Fvec,-1,varargin{:});  % compute JAC'*Fvec
normgradinf = norm(grad,inf);

% Print header.
header = sprintf(['\n                                         Norm of      First-order   Trust-region\n',...
                    ' Iteration  Func-count     f(x)          step         optimality    radius']);
formatstr = ' %5.0f      %5.0f   %13.6g  %13.6g   %12.3g    %12.3g';
if verbosity > 1
  disp(header);
end

% Initialize the output function.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = systemSolve.systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'init',iter, ...
        numFevals,Fvec,obj,[],[],[],Delta,stepAccept,varargin{:});
    if stop
        [x,Fvec,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.systemSolve.cleanUpInterrupt(xOutputfcn,optimValues);
        return;
    end
end

% If using Broyden updates need to provide an initial inverse Jacobian.
if broyden, 
  ws = warning('off');
  invJAC = JAC\speye(nfnc);
  warning(ws);
else
  invJAC = [];
end

% Compute initial diagonal scaling matrix.
if scale
  if giventypx && ~isempty(typx) % scale based on typx values
    typx(typx==0) = 1; % replace any zero entries with ones
    scalMat = 1./abs(typx);
  else         % scale based on norm of the Jacobian (not currently active)  
    scalMat = systemSolve.getscalMat(nvar,JAC,scalemin,scalemax);
  end
end

% Display initial iteration information.
formatstr0 = ' %5.0f      %5.0f   %13.6g                  %12.3g    %12.3g';
% obj is 0.5*F'*F but want to display F'*F
iterOutput0 = sprintf(formatstr0,iter,numFevals,2*obj,normgradinf,Delta);
if verbosity > 1
   disp(iterOutput0);
end
% OutputFcn call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = systemSolve.systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'iter',iter, ...
        numFevals,Fvec,obj,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
    if stop
        [x,Fvec,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.systemSolve.cleanUpInterrupt(xOutputfcn,optimValues);
        return;
    end
end


% Test convergence at initial point.
[exitStatus,done,EXITFLAG,msg] = systemSolve.gtestStop(normgradinf,tolf,tolx,...
     stepAccept,iter,maxit,numFevals,maxfunc,Delta,normd,...
     obj,objold,d,xvec,broydenStep);

% Beginning of main iteration loop.
while ~done
  iter = iter + 1;

  % Compute step, d, using dogleg approach.
  [d,quadObj,normd,normdscal,illcondition] = ...
       systemSolve.dogleg(nvar,nfnc,Fvec,JAC,grad,Delta,d,invJAC,broyden, ...
       scalMat,mtxmpy,varargin);
  if broydenJac 
    broydenStep = true; 
  else
    broydenStep = false;
  end

  % Compute the model reduction given by d (pred).
  pred = -quadObj;

  % Compute the trial point, xTrial.
  xTrial = xvec + d;

  % Evaluate nonlinear equations and objective at trial point.
  x(:) = xTrial; % reshape xTrial to a matrix for evaluations. 
  switch funfcn{1}
  case 'fun'
    F = feval(funfcn{3},x,varargin{:});
  case 'fungrad'
    [F,JACTrial] = feval(funfcn{3},x,varargin{:});
    numJevals = numJevals + 1;
  case 'fun_then_grad'
    F = feval(funfcn{3},x,varargin{:}); 
  otherwise
    error('optim:trustnleqn:UndefinedCalltype','Undefined calltype in FSOLVE.')
  end  
  numFevals = numFevals + 1;
  FTrial = F(:); % make FTrial a vector
  objTrial = 0.5*FTrial'*FTrial; 

  % Compute the actual reduction given by xTrial (ared).
  ared = obj - objTrial;

  % Compute ratio = ared/pred.
  if pred <= 0 % reject step
    ratio = 0;
  else
    ratio = ared/pred;
  end
  
  % OutputFcn call
  if haveoutputfcn
      stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,'interrupt',varargin{:});
      if stop  % Stop per user request.
          [x,Fvec,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.systemSolve.cleanUpInterrupt(xOutputfcn,optimValues);
          return;
      end
  end
  
  if ratio > eta1 % accept step.

    % Update information.
    if broyden, numReject = 0; Fdiff = FTrial - Fvec; end
    xvec = xTrial; Fvec = FTrial; objold = obj; obj = objTrial;
    x(:) = xvec; % update matrix representation

    % Compute JAC at new point. (already computed with F if 'fungrad')

    % Broyden update.
    if broyden && ~illcondition
      normd2 = normd^2;
      denom = d'*(invJAC*Fdiff);
      if abs(denom) < 0.1*normd2
        alpha = 0.8;
        denom = alpha*denom + (1-alpha)*normd2;
      else
        alpha = 1.0;
      end  
      JACd = feval(mtxmpy,JAC,d,1,varargin{:});  % compute JAC*d        
      JACbroyden = JAC + alpha*((Fdiff - JACd)*d')/(normd^2);
      invJAC = invJAC + alpha*((d - invJAC*Fdiff)*(d'*invJAC))/denom;
    end

    % Compute sparse finite difference Jacobian if needed.
    if ~gradflag && (~broyden || (broyden && illcondition))
      [JACfindiff,numFDfevals] = systemSolve.sfdnls(x,Fvec,JACfindiff,group,[], ...
                  DiffMinChange,DiffMaxChange,funfcn{3},[],[],varargin{:});
      numFevals = numFevals + numFDfevals;
    end

    if broyden && ~illcondition
      JAC = JACbroyden; broydenJac = true;
    else
      switch funfcn{1}
      case 'fun'
        JAC = JACfindiff;
      case 'fungrad'
        JAC = JACTrial;
      case 'fun_then_grad'
        JAC = feval(funfcn{4},x,varargin{:});
        numJevals = numJevals + 1;
      otherwise
        error('optim:trustnleqn:UndefinedCalltype','Undefined calltype in FSOLVE.')
      end
      broydenJac = false;
    end 
    grad = feval(mtxmpy,JAC,Fvec,-1,varargin{:});  % compute JAC'*Fvec     
    normgradinf = norm(grad,inf);

    if broyden && illcondition % update inverse Jacobian
      ws = warning('off');
      invJAC = JAC\speye(nfnc);
      warning(ws);
    end

    % Update internal diagonal scaling matrix (dynamic scaling).
    if scale && ~giventypx
      scalMat = systemSolve.getscalMat(nvar,JAC,scalemin,scalemax);
    end

    stepAccept = true;

  else % reject step.
 
    if broyden, numReject = numReject + 1; end
    stepAccept = false;

  end 

  % Print iteration statistics.
  if verbosity > 1
      % obj is 0.5*F'*F but want to display F'*F
      iterOutput = sprintf(formatstr,iter,numFevals,2*obj,normd,normgradinf,Delta);
      disp(iterOutput);
  end
  % OutputFcn call
  if haveoutputfcn || haveplotfcn
      [xOutputfcn, optimValues, stop] = systemSolve.systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'iter',iter, ...
        numFevals,Fvec,obj,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
      if stop
          [x,Fvec,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.systemSolve.cleanUpInterrupt(xOutputfcn,optimValues);
          return;
      end
  end

  % Update trust region radius.
  Delta = systemSolve.updateDelta(Delta,ratio,normdscal,eta1,eta2,...
                      alpha1,alpha2,DeltaMax);

  % Check for termination.
  [exitStatus,done,EXITFLAG,msg] = systemSolve.gtestStop(normgradinf,tolf,tolx,...
       stepAccept,iter,maxit,numFevals,maxfunc,Delta,normd,...
       obj,objold,d,xvec,broydenStep);
  
  % If using Broyden updating to compute the Jacobian and the last two steps
  % were rejected or the Broyden Jacobian is ill-conditioned this may indicate 
  % that the Jacobian approximation is not accurate.  Compute a new finite 
  % difference or analytic Jacobian.  
  if broydenJac
    if (numReject == 2 && ~done) || (illcondition && ~done)
      EXITFLAG = 0;
      % Compute new finite difference or analytic Jacobian.
      x(:) = xvec; % reset matrix representation of x to current value in xvec
      if ~gradflag
        [JACfindiff,numFDfevals] = sfdnls(x,Fvec,JACfindiff,group,[], ...
                   DiffMinChange,DiffMaxChange,funfcn{3},[],[],varargin{:});
        numFevals = numFevals + numFDfevals;
      end
      switch funfcn{1}
      case 'fun'
        JAC = JACfindiff;
      case 'fungrad'
        [F,JAC] = feval(funfcn{3},x,varargin{:});
        numFevals = numFevals + 1;
        numJevals = numJevals + 1;
      case 'fun_then_grad'
        JAC = feval(funfcn{4},x,varargin{:});
        numJevals = numJevals + 1;
      otherwise
        error('optim:trustnleqn:UndefinedCalltype','Undefined calltype in FSOLVE.')
      end
      grad = feval(mtxmpy,JAC,Fvec,-1,varargin{:});  % compute JAC'*Fvec
      normgradinf = norm(grad,inf);
      broydenJac = false;
      ws = warning('off');
      invJAC = JAC\speye(nfnc);
      warning(ws);
      if scale && ~giventypx     % Update internal scaling matrix.
        scalMat = systemSolve.getscalMat(nvar,JAC,scalemin,scalemax);
      end
    end
  end
end

if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues] = systemSolve.systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'done',iter, ...
        numFevals,Fvec,obj,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
    % Optimization done, so ignore "stop"
end


% Optimization is finished.

% Assign output statistics.
OUTPUT.iterations = iter;
OUTPUT.funcCount = numFevals;
OUTPUT.algorithm = 'trust-region dogleg';
OUTPUT.firstorderopt = normgradinf;

% TRUSTNLEQN finished
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxfunc,maxit,tolf,tolx,derivCheck,DiffMinChange,...
          DiffMaxChange,mtxmpy,typx,giventypx,JACfindiff,structure,outputfcn,plotfcns] = ...
          getOpts(nfnc,nvar,options,defaultopt,gradflag)
%getOpts gets the user-defined options for TRUSTNLEQN.

% Both Medium and Large-Scale options.
maxfunc = optimget(options,'MaxFunEvals',defaultopt,'fast');
if ischar(maxfunc)
  if isequal(lower(maxfunc),'100*numberofvariables')
    maxfunc = 100*nvar;
  else
    error('optim:trustnleqn:InvalidMaxFunEvals', ...
          'Option ''MaxFunEvals'' must be an integer value if not the default.')
  end
end
maxit = optimget(options,'MaxIter',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
if isfield(options,'OutputFcn')
    outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
else
    outputfcn = defaultopt.OutputFcn;
end

if isfield(options,'PlotFcns')
    plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
else
    plotfcns = defaultopt.PlotFcns;
end

% Medium-Scale only options.
derivCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');

% Large-Scale only options.
mtxmpy = optimget(options,'JacobMult',defaultopt,'fast');
if isempty(mtxmpy)
  mtxmpy = @systemSolve.atamult;
end
giventypx = true;
typx = optimget(options,'TypicalX',defaultopt,'fast');
if ischar(typx)
  if isequal(lower(typx),'ones(numberofvariables,1)')
    typx = ones(nvar,1);
    giventypx = false;
  else
    error('optim:trustnleqn:InvalidTypicalX', ...
          'Option ''TypicalX'' must be a matrix (not a string) if not the default.')
  end
end
structure = true;
if ~gradflag
  JACfindiff = optimget(options,'JacobPattern',defaultopt,'fast');
  if ischar(JACfindiff) 
    if isequal(lower(JACfindiff),'sparse(ones(jrows,jcols))')
      JACfindiff = sparse(ones(nfnc,nvar));
      structure = false;
    else
      error('optim:trustnleqn:InvalidJacobPattern', ...
            'Option ''JacobPattern'' must be a matrix if not the default.')
    end
  end
else
  JACfindiff = [];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [exitStatus,done,EXITFLAG,msg] = gtestStop(normgradinf,tolf,tolx,...
     stepAccept,iter,maxit,numFevals,maxfunc,Delta,normd,...
     obj,objold,d,xvec,broydenStep)
%gtestStop checks the termination criteria for TRUSTNLEQN.

exitStatus = 0;
done = false;
EXITFLAG = 0;
msg = '';

% Check termination criteria.
if stepAccept && normgradinf < tolf
  msg = sprintf('Optimization terminated: first-order optimality is less than options.TolFun.');
  exitStatus = 1;
  done = true;
  EXITFLAG = 1;
elseif iter > 1 && max(abs(d)./(abs(xvec)+1)) < max(tolx^2,eps) 
   if 2*obj < sqrt(tolf) % fval'*fval < sqrt(tolf)
      msg = sprintf(['Optimization terminated: norm of relative change in X is less\n' ...
                     ' than max(options.TolX^2,eps) and  sum-of-squares of function \n' ...
                     ' values is less than sqrt(options.TolFun).']);
      exitStatus = 2;
      EXITFLAG = 2;
   else
      msg = sprintf(['Optimizer appears to be converging to a point which is not a root.\n',...
         ' Norm of relative change in X is less than max(options.TolX^2,eps) but\n',...
         ' sum-of-squares of function values is greater than or equal to sqrt(options.TolFun)\n',...  
         ' Try again with a new starting guess.']);
      exitStatus = 3;
      EXITFLAG = -2;
    end
    done = true;
elseif iter > 1 && stepAccept && normd < 0.9*Delta ...
                && abs(objold-obj) < max(tolf^2,eps)*(1+abs(objold))
  if 2*obj < sqrt(tolf) % fval'*fval < sqrt(tolf)
     msg = sprintf(['Optimization terminated: relative function value changing by less\n' ...
                    ' than max(options.TolFun^2,eps) and sum-of-squares of function\n' ... 
                    ' values is less than sqrt(options.TolFun).']);
     exitStatus = 4;
     EXITFLAG = 3;
  else
      msg = sprintf(['Optimizer appears to be converging to a point which is not a root.\n',...
         ' Relative function value changing by less than max(options.TolFun^2,eps) but\n',...
         ' sum-of-squares of function values is greater than or equal to sqrt(options.TolFun)\n',...  
        ' Try again with a new starting guess.']);
      exitStatus = 5;
      EXITFLAG = -2; 
  end
  done = true;
elseif Delta < 2*eps
  msg = sprintf(['Optimization terminated: no further progress can be made.\n',...
         ' Trust-region radius less than 2*eps.\n',...
         ' Problem may be ill-conditioned or Jacobian may be inaccurate.\n',...
         ' Try using exact Jacobian or check Jacobian for errors.']);
  exitStatus = 6;
  done = true;
  EXITFLAG = -3;
elseif iter >= maxit
  msg = sprintf(['Maximum number of iterations reached:\n',...
                 ' increase options.MaxIter.']);
  exitStatus = 7;
  done = true;
  EXITFLAG = 0;
elseif numFevals >= maxfunc
  msg = sprintf(['Maximum number of function evaluations reached:\n',...
                 ' increase options.MaxFunEvals.']);
  exitStatus = 8;
  done = true;
  EXITFLAG = 0;  
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Delta = updateDelta(Delta,ratio,normdscal,eta1,eta2,...
                             alpha1,alpha2,DeltaMax)
%updateDelta updates the trust region radius in TRUSTNLEQN.
%
%   systemSolve.updateDelta updates the trust region radius based on the value of
%   ratio and the norm of the scaled step.

if ratio < eta1
  Delta = alpha2*normdscal;
elseif ratio >= eta2
  Delta = max(Delta,alpha1*normdscal);
end
Delta = min(Delta,DeltaMax);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scalMat = getscalMat(nvar,JAC,scalemin,scalemax)
%systemSolve.getscalMat computes the scaling matrix in TRUSTNLEQN.
%
%   getscalMat computes the scaling matrix based on the norms 
%   of the columns of the Jacobian.

scalMat = ones(nvar,1);
for i=1:nvar
  scalMat(i,1) = norm(JAC(:,i));
end
scalMat(scalMat<scalemin) = scalemin;  % replace small entries
scalMat(scalMat>scalemax) = scalemax;  % replace large entries
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,state,iter,numFevals, ...
    Fvec,obj,normd,grad,normgradinf,Delta,stepAccept,varargin)
% systemSolve.callOutputAndPlotFcns assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter', or 'done'. 
% We do not handle the case 'interrupt' because we do not want to update
% xOutputfcn or optimValues (since the values could be inconsistent) before calling
% the outputfcn; in that case the outputfcn is called directly rather than
% calling it inside systemSolve.systemSolve.callOutputAndPlotFcns.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.

optimValues.iteration = iter;
optimValues.funccount = numFevals;
optimValues.fval = Fvec;
optimValues.stepsize = normd; 
optimValues.gradient = grad; 
optimValues.firstorderopt = normgradinf;
optimValues.trustregionradius = Delta;
optimValues.stepaccept = stepAccept;

xOutputfcn(:) = xvec;  % Set xvec to have user expected size
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('optim:trustnleqn:UnknownStateInsystemSolve.systemSolve.callOutputAndPlotFcns','Unknown state in systemSolve.systemSolve.callOutputAndPlotFcns.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('optim:trustnleqn:UnknownStateInsystemSolve.systemSolve.callOutputAndPlotFcns','Unknown state in systemSolve.systemSolve.callOutputAndPlotFcns.')
    end
end
end
%--------------------------------------------------------------------------
function [x,Fvec,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.systemSolve.cleanUpInterrupt(xOutputfcn,optimValues)
% systemSolve.systemSolve.cleanUpInterrupt sets the outputs arguments to be the values at the last call
% of the outputfcn during an 'iter' call (when these values were last known to
% be consistent). 

x = xOutputfcn; 
Fvec = optimValues.fval;
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'trust-region dogleg';
OUTPUT.firstorderopt = optimValues.firstorderopt; 
JAC = []; % May be in an inconsistent state
msg = 'Optimization terminated prematurely by user.';


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J,ncol] = sfdnls(xcurr,valx,Jstr,group,alpha,DiffMinChange, ...
                            DiffMaxChange,fun,XDATA,YDATA,varargin)
%SFDNLS    Sparse Jacobian via finite differences
%
% J = sfdnls(xcurr,valx,Jstr,group,[],DiffMinChange,DiffMaxChange,fun, ...
% YDATA,varargin) returns the sparse finite difference approximation J of 
% the Jacobian matrix of the function 'fun' at the current point xcurr. The 
% vector group indicates how to use sparse finite differencing: group(i) = j 
% means that column i belongs to group (or color) j. Each group (or color) 
% corresponds to a function difference. The input varargin contains the extra 
% parameters (possibly) needed by function 'fun'. 
% DiffMinChange and DiffMaxChange indicate, respectively, the minimum and 
% maximum change in variables during the finite difference calculation.
%
% A non-empty input alpha overrides the default finite differencing stepsize.
%
% [J,ncol] = sfdnls(...) returns the number of function evaluations used
% in ncol.

%   Copyright 1990-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:43 $

%
if nargin < 8
   error('optim:sfdnls:RequiresEightArguments','SFDNLS requires at least eight arguments.')
elseif nargin < 10
   YDATA = [];
elseif nargin < 9
   XDATA = [];
end

scalealpha = false;
x = xcurr(:); % make it a vector
[m,n] = size(Jstr); 
ncol = max(group); 
if isempty(alpha)
    scalealpha = true;
    alpha = repmat(sqrt(eps),ncol,1);
end
J = spones(Jstr);

% If lsqcurvefit, then add XDATA to objective's input list.
% xargin{1} will be updated right before each evaluation
if ~isempty(XDATA)
    xargin = {xcurr,XDATA};
else
    xargin = {xcurr};
end

if ncol < n
   for k = 1:ncol
      d = (group == k);
      if scalealpha
         xnrm = norm(x(d));
         xnrm = max(xnrm,1);
         alpha(k) = alpha(k)*xnrm;
      end
      
      % Ensure magnitude of step-size lies within interval 
      % [DiffMinChange, DiffMaxChange]
      alpha(k) = sign(alpha(k))*min(max(abs(alpha(k)),DiffMinChange), ...
                                  DiffMaxChange);      
      y = x + alpha(k)*d;
      
      xcurr(:) = y;  % reshape for userfunction
      xargin{1} = xcurr; % update x in list of input arguments to objective
      v = feval(fun,xargin{:},varargin{:});
      if ~isempty(YDATA)
         v = v - YDATA;
      end
      v = v(:);
      
      w = (v-valx)/alpha(k);
      cols = find(d); 
      
      A = sparse(m,n);
      A(:,cols) = J(:,cols);
      J(:,cols) = J(:,cols) - A(:,cols);
      [i,j,val] = find(A);
      [p,ind] = sort(i);
      val(ind) = w(p);
      A = sparse(i,j,full(val),m,n);
      J = J + A;
   end
else % ncol ==n
   J = full(J);
   for k = 1:n
      if scalealpha
         xnrm = norm(x(k));
         xnrm = max(xnrm,1);
         alpha(k) = alpha(k)*xnrm;
      end
      
      % Ensure magnitude of step-size lies within interval 
      % [DiffMinChange, DiffMaxChange]
      alpha(k) = sign(alpha(k))*min(max(abs(alpha(k)),DiffMinChange), ...
                                  DiffMaxChange);      
      y = x;
      y(k) = y(k) + alpha(k);

      xcurr(:) = y;  % reshape for userfunction
      xargin{1} = xcurr; % update x in list of input arguments to objective
      v = feval(fun,xargin{:},varargin{:});
      if ~isempty(YDATA)
         v = v - YDATA;
      end
      v = v(:);
      J(:,k) = (v-valx)/alpha(k);
   end
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[V] = atamult(A,Y,flag,varargin)
%ATAMULT Jacobian-matrix multiply
%
%	V = ATAMULT(A,Y) computes V = (A'*(A*Y)).
%
%	V = ATAMULT(A,Y,flag) computes V = (A'*(A*Y)) if flag = 0,
%                                  V = A*Y        if flag > 0,
%                                  V = A'*Y       if flag < 0.
%
% Note: varargin is not used but must be provided in case 
% the objective function has additional problem dependent
% parameters (which will be passed to this routine as well).

%   Copyright 1990-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:47:15 $

if nargin < 3 || flag == 0
   V = (A'*(A*Y));
elseif flag > 0
   V = A*Y;
else
   V = A'*Y;
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,quadObj,normd,normdscal,illcondition] = ...
     dogleg(nvar,nfnc,F,JAC,grad,Delta,d,invJAC,broyden, ...
     scalMat,mtxmpy,varargin)
%DOGLEG approximately solves trust region subproblem via a dogleg approach.
%
%   DOGLEG finds an approximate solution d to the problem:
%
%     min_d      f + g'd + 0.5*d'Bd
%
%     subject to ||Dd|| <= Delta
%
%   where g is the gradient of f, B is a Hessian approximation of f, D is a 
%   diagonal scaling matrix and Delta is a given trust region radius.    

%   Copyright 1990-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2006/12/15 19:28:59 $
%   Richard Waltz, June 2001
%
% NOTE: The scaling matrix D above is called scalMat in this routine.

illcondition = 0;

% Initialize local arrays.
dCauchy   = zeros(nvar,1);
dGN       = zeros(nvar,1);
gradscal  = zeros(nvar,1);
gradscal2 = zeros(nvar,1);
JACvec    = zeros(nfnc,1);

% Compute scaled gradient and other scaled terms.
gradscal = grad./scalMat;
gradscal2 = gradscal./scalMat;
normgradscal = norm(gradscal);

if normgradscal >= eps

  % First compute the Cauchy step (in scaled space).
  dCauchy = -(Delta/normgradscal)*gradscal;
  JACvec = feval(mtxmpy,JAC,gradscal2,1,varargin{:});
  denom = Delta*JACvec'*JACvec;
  tauterm = normgradscal^3/denom;
  tauC = min(1,tauterm);
  dCauchy = tauC*dCauchy;

  % Compute quadratic objective at Cauchy point.
  JACvec = feval(mtxmpy,JAC,dCauchy./scalMat,1,varargin{:});  % compute JAC*d
  objCauchy = gradscal'*dCauchy + 0.5*JACvec'*JACvec;

  normdCauchy = min(norm(dCauchy),Delta);

else

  % Set Cauchy step to zero step and continue.
  objCauchy = 0;
  normdCauchy = 0;

end  

if Delta - normdCauchy < eps;

  % Take the Cauchy step if it is at the boundary of the trust region.
  d = dCauchy; quadObj = objCauchy;

else

  condition = condest(JAC);
  if condition > 1.0e10, illcondition = 1; end

  if condition > 1.0e15

    % Take the Cauchy step if Jacobian is (nearly) singular. 
    d = dCauchy; quadObj = objCauchy;

  else

    % Compute the Gauss-Newton step (in scaled space).
    if broyden
      dGN = -invJAC*F;
    else
      % Disable the warnings about conditioning for singular and
      % nearly singular matrices
      warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
      warningstate2 = warning('off', 'MATLAB:singularMatrix');
      dGN = -JAC\F;
      % Restore the warning states to their original settings
      warning(warningstate1)
      warning(warningstate2)
    end
    dGN = dGN.*scalMat;     % scale the step
    
    if any(~isfinite(dGN))  

      % Take the Cauchy step if the Gauss-Newton step gives bad values.
      d = dCauchy; quadObj = objCauchy;
    else

      normdGN = norm(dGN);
      if normdGN <= Delta

        % Compute quadratic objective at Gauss-Newton point.
        JACvec = feval(mtxmpy,JAC,dGN./scalMat,1,varargin{:});  % compute JAC*d
        objGN = gradscal'*dGN + 0.5*JACvec'*JACvec;

        if ~illcondition
          % Take Gauss-Newton step if inside trust region and well-conditioned.
          d = dGN; quadObj = objGN; 
        else
          % Compare Cauchy step and Gauss-Newton step if ill-conditioned.
          if objCauchy < objGN
            d = dCauchy; quadObj = objCauchy;
          else
            d = dGN; quadObj = objGN;
          end
        end        

      else

        % Find the intersect point along dogleg path.

        Delta2 = Delta^2;
        normdCauchy2 = min(normdCauchy^2,Delta2);
        normdGN2 = normdGN^2;
        dCdGN = dCauchy'*dGN;
        dCdGNdist2 = max((normdCauchy2+normdGN2-2*dCdGN),0);

        if dCdGNdist2 == 0
          tauI = 0;
        else 
          tauI = (normdCauchy2-dCdGN + sqrt((dCdGN-normdCauchy2)^2 ...
                + dCdGNdist2*(Delta2-normdCauchy2))) / dCdGNdist2;
        end

        d = dCauchy + tauI*(dGN-dCauchy);

        % Compute quadratic objective at intersect point.
        JACvec = feval(mtxmpy,JAC,d./scalMat,1,varargin{:});  % compute JAC*d
        objIntersect = gradscal'*d + 0.5*JACvec'*JACvec;        

        if ~illcondition
          % Take Intersect step if well-conditioned.
          quadObj = objIntersect; 
        else
          % Compare Cauchy step and Intersect step if ill-conditioned.
          if objCauchy < objIntersect
            d = dCauchy; quadObj = objCauchy;
          else
            quadObj = objIntersect;
          end
        end            

      end
    end
  end    
end

% The step computed was the scaled step.  Unscale it.
normdscal = norm(d);
d = d./scalMat;
normd = norm(d);




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function step = cubici2(c,f,x)
%CUBICI2 Determine optimizer step from three points and one gradient.
%   STEP = CUBICI2(c,f,x)
%   Finds the cubic p(x) with p(x(1:3)) = f(1:3) and p'(0) = c.
%   Returns the minimizer of p(x) if it is positive.
%   Calls QUADI if the minimizer is negative.

%   Copyright 1990-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:11 $

% p(x) = a/3*x^3 - b*x^2 + c*x + d.
% c = p'(0) is the first input parameter.
% Solve [1/3*x.^3 -1/2*x^2 ones(3,1)]*[a b d]' = f - c*x.
% Compute a and b; don't need d.
%    a = 3*(x1^2*(f2-f3) + x2^2*(f3-f1) + x3^2*(f1-f2))/h
%    b = (x1^3*(f2-f3) + x2^3*(f3-f1) + x3^3*(f1-f2))/h
%    where h = (x1-x2)*(x2-x3)*(x3-x1)*(x1*x2 + x2*x3 + x3*x1).
% Local min and max where p'(s) = a*s^2 - 2*b*s + c = 0
% Local min always comes from plus sign in the quadratic formula.
% If p'(x) has no real roots, step = b/a.
% If step < 0, use quadi instead.

x = x(:);
f = f(:);
g = f - c*x;
g = g([2 3 1]) - g([3 1 2]);
y = x([2 3 1]);
h = prod(x-y)*(x'*y);
a = 3*(x.^2)'*g/h;
b = (x.^3)'*g/h;

% Find minimizer.
step = (b + real(sqrt(b^2-a*c)))/a;

% Is step acceptable?
if step < 0 | ~isfinite(step)
   step = abs(systemSolve.quadi(x,f));
end
if isnan(step)
   step = x(2)/2;
end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,f] = cubici3(f2,f1,c2,c1,dx)
%CUBICI3  Cubicly interpolates 2 points and gradients to find step and min.
%   This function uses cubic interpolation and the values of 
%   two points and their gradients in order to estimate the minimum s of a 
%   a function along a line and returns s and f=F(s);
%

%  The equation is F(s) = a/3*s^3 + b*s^2 + c1*s + f1
%      and F'(s) = a*s^2+2*b*s + c1
%  where we know that 
%          F(0) = f1
%          F'(0) = c1  
%          F(dx) = f2   implies: a/3*dx^3 + b*dx^2 + c1*dx + f1 = f2
%          F'(dx) = c2  implies: a*dx^2+2*b*dx + c1 = c2

%   Copyright 1990-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:12 $

if isinf(f2), 
    f2 = 1/eps; 
end
a = (6*(f1-f2)+3*(c1+c2)*dx)/dx^3;
b = (3*(f2-f1)-(2*c1+c2)*dx)/dx^2;
disc = b^2 - a*c1;
if a==0 & b==0 
    % F(s) is linear: F'(s) = c1, which is never zero;
    % minimum is s=Inf or -Inf (we return s always positive so s=Inf).
    s = inf; 
elseif a == 0 
    % F(s) is quadratic so we know minimum s
    s = -c1/(2*b);
elseif disc <= 0
    % If disc = 0 this is exact. 
    % If disc < 0 we ignore the complex component of the root.
    s = -b/a;  
else
    s = (-b+sqrt(disc))/a;
end
if s<0,  s = -s; end
if isinf(s)
    f = inf;
else
    % User Horner's rule
    f = ((a/3*s + b)*s + c1)*s + f1;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s] = finitedifferences(xCurrent,...
             xOriginalShape,funfcn,confcn,lb,ub,fCurrent,cCurrent,XDATA,YDATA,...
             DiffMinChange,DiffMaxChange,typicalx,finDiffType,variables,...
             LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,isFseminf,varargin)
%FINITEDIFFERENCES computes finite-difference derivatives.
%
% This helper function computes finite-difference derivatives of the objective 
% and constraint functions.
%
%  [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s] = FINITEDIFFERENCES(xCurrent, ... 
%                  xOriginalShape,funfcn,confcn,lb,ub,fCurrent,cCurrent, ...
%                  XDATA,YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType, ...
%                  variables,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s, ...
%                  varargin)
% computes the finite-difference gradients of the objective and
% constraint functions.
%
%  gradf = FINITEDIFFERENCES(xCurrent,xOriginalShape,funfcn,[],lb,ub,fCurrent,...
%                  [],YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
%                  variables,[],[],[],[],[],[],varargin)
% computes the finite-difference gradients of the objective function.
%
%
% INPUT:
% xCurrent              Point where gradient is desired
% xOriginalShape        Shape of the vector of variables supplied by the user
%                       (The value of xOriginalShape is NOT used)
% funfcn, confcn        Cell arrays containing info about objective and
%                       constraints, respectively. The objective (constraint) 
%                       derivatives are computed if and only if funfcn 
%                       (confcn) is nonempty.
%                       
% lb, ub                Lower and upper bounds
% fCurrent, cCurrent    Values at xCurrent of the function and the constraints 
%                       to be differentiated. Note that fCurrent can be a scalar 
%                       or a (row or column) vector. 
%
% XDATA, YDATA          Data passed from lsqcurvefit.
%
% DiffMinChange, 
% DiffMaxChange         Minimum and maximum values of perturbation of xCurrent 
% finDiffType           Type of finite difference desired (only forward 
%                       differences implemented so far)
% variables             Variables w.r.t which we want to differentiate. Possible 
%                       values are 'all' or an integer between 1 and the
%                       number of variables.
%
% LAMBDA,NEWLAMBDA,
% OLDLAMBDA,POINT,
% FLAG,s                Parameters for semi-infinite constraints
%
% isFseminf             True if caller is fseminf, false otherwise
%
% varargin              Problem-dependent parameters passed to the objective and 
%                       constraint functions
%
% OUTPUT:
% gradf                 If fCurrent is a scalar, gradf is the finite-difference 
%                       gradient of the objective; if fCurrent is a vector,
%                       gradf is the finite-difference Jacobian  
% cJac                  Finite-difference Jacobian of the constraints
% NEWLAMBDA,
% OLDLAMBDA,s           Parameters for semi-infinite constraints

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.4.8 $  $Date: 2006/12/15 19:29:32 $

% For vector-valued functions in funfcn, we make the function values 
% (both at xCurrent and those at the perturbed points) column vectors 
% to ensure that the given fCurrent and the computed fplus will have the 
% same shape so that fplus - fCurrent be well defined.  
fCurrent = fCurrent(:);
numberOfFunctions = numel(fCurrent);
numberOfVariables = numel(xCurrent); 
functionIsScalar = (numberOfFunctions == 1);

% nonEmptyLowerBounds = true if lb is not empty, false if it's empty;
% analogoulsy for nonEmptyUpperBound
nonEmptyLowerBounds = ~isempty(lb);
nonEmptyUpperBounds = ~isempty(ub);

% Make sure xCurrent and typicalx are column vectors so that the 
% operation max(abs(xCurrent),abs(typicalx)) won't error
xCurrent = xCurrent(:); typicalx = typicalx(:);
% Value of stepsize suggested in Trust Region Methods, Conn-Gould-Toint, section 8.4.3
CHG = sqrt(eps)*sign(xCurrent).*max(abs(xCurrent),abs(typicalx));
%
% Make sure step size lies within DiffminChange and DiffMaxChange
%
CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
len_cCurrent = length(cCurrent); % For semi-infinite

if nargout < 3
   NEWLAMBDA=[]; OLDLAMBDA=[]; s=[];
end
if nargout > 1
      cJac = zeros(len_cCurrent,numberOfVariables);  
else
      cJac = [];
end
% allVariables = true/false if finite-differencing wrt to all/one variables
allVariables = false;
if ischar(variables)
   if strcmp(variables,'all')
      variables = 1:numberOfVariables;
      allVariables = true;
   else
      error('optimlib:finitedifferences:InvalidVariables', ...
            'Unknown value of input ''variables''.')
   end
end

% Preallocate gradf for speed 
if ~isempty(funfcn)
    if functionIsScalar
        gradf = zeros(numberOfVariables,1);
    elseif allVariables % vector-function and gradf estimates full Jacobian
        gradf = zeros(numberOfFunctions,numberOfVariables);
    else % vector-function and gradf estimates one column of Jacobian
        gradf = zeros(numberOfFunctions,1);
    end
else
    gradf = [];
end

% Do this switch outside of loop for speed
if isFseminf
   vararginAdjusted = varargin(3:end);
else
   vararginAdjusted = varargin;     
end

% If lsqcurvefit, then add XDATA to objective's input list.
% xargin{1} will be updated right before each evaluation 
if ~isempty(XDATA)
    xargin = {xOriginalShape,XDATA};
else
    xargin = {xOriginalShape};
end

for gcnt=variables
   if gcnt == numberOfVariables, 
      FLAG = -1; 
   end
   temp = xCurrent(gcnt);
   xCurrent(gcnt)= temp + CHG(gcnt);
         
   if (nonEmptyLowerBounds && isfinite(lb(gcnt))) || (nonEmptyUpperBounds && isfinite(ub(gcnt)))
      % Enforce bounds while finite-differencing.
      % Need lb(gcnt) ~= ub(gcnt), and lb(gcnt) <= temp <= ub(gcnt) to enforce bounds.
      % (If the last qpsub problem was 'infeasible', the bounds could be currently violated.)
      if (lb(gcnt) ~= ub(gcnt)) && (temp >= lb(gcnt)) && (temp <= ub(gcnt)) 
          if  ((xCurrent(gcnt) > ub(gcnt)) || (xCurrent(gcnt) < lb(gcnt))) % outside bound ?
              CHG(gcnt) = -CHG(gcnt);
              xCurrent(gcnt)= temp + CHG(gcnt);
              if (xCurrent(gcnt) > ub(gcnt)) || (xCurrent(gcnt) < lb(gcnt)) % outside other bound ?
                  [newchg,indsign] = max([temp-lb(gcnt), ub(gcnt)-temp]);  % largest distance to bound
                  if newchg >= DiffMinChange
                      CHG(gcnt) = ((-1)^indsign)*newchg;  % make sure sign is correct
                      xCurrent(gcnt)= temp + CHG(gcnt);

                      % This warning should only be active if the user doesn't supply gradients;
                      % it shouldn't be active if we're doing derivative check 
                      warning('optimlib:finitedifferences:StepReduced', ...
                             ['Derivative finite-differencing step was artificially reduced to be within\n', ...
                              'bound constraints. This may adversely affect convergence. Increasing distance between\n', ...
                              'bound constraints, in dimension %d to be at least %0.5g may improve results.'], ...
                              gcnt,abs(2*CHG(gcnt)))
                  else
                      error('optimlib:finitedifferences:DistanceTooSmall', ...
                          ['Distance between lower and upper bounds, in dimension %d is too small to compute\n', ...
                           'finite-difference approximation of derivative. Increase distance between these\n', ...
                           'bounds to be at least %0.5g.'],gcnt,2*DiffMinChange)
                  end          
              end
          end
      end
   end % of 'if isfinite(lb(gcnt)) || isfinite(ub(gcnt))'
   
   xOriginalShape(:) = xCurrent;
   xargin{1} = xOriginalShape; % update x in list of input arguments to objective
   if ~isempty(funfcn) % Objective gradient required
       % The length of varargin (depending on the caller being fseminf or not)
       % was calculated outside of the loop for speed.
       fplus = feval(funfcn{3},xargin{:},vararginAdjusted{:});
   else
       fplus = [];
   end
   % YDATA: Only used by lsqcurvefit, which has no nonlinear constraints
   % (the only type of constraints we do finite differences on: bounds 
   % and linear constraints do not require finite differences) and thus 
   % no needed after evaluation of constraints
   if ~isempty(YDATA)    
      fplus = fplus - YDATA;
   end
   % Make sure it's in column form
   fplus = fplus(:);
   if ~isempty(funfcn)
       if functionIsScalar
           gradf(gcnt,1) =  (fplus-fCurrent)/CHG(gcnt);
       elseif allVariables % vector-function and gradf estimates full Jacobian
           gradf(:,gcnt) = (fplus-fCurrent)/CHG(gcnt);
       else % vector-function and gradf estimates only one column of Jacobian
           gradf = (fplus-fCurrent)/CHG(gcnt);
       end
   end
   
   if ~isempty(cJac) % Constraint gradient required
         if isFseminf 
            [ctmp,ceqtmp,NPOINT,NEWLAMBDA,OLDLAMBDA,LOLD,s] = ...
               semicon(xOriginalShape,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});            
         else
            [ctmp,ceqtmp] = feval(confcn{3},xOriginalShape,varargin{:});
         end
         cplus = [ceqtmp(:); ctmp(:)];

      % Next line used for problems with varying number of constraints
      if isFseminf && len_cCurrent~=length(cplus)
         cplus=v2sort(cCurrent,cplus); 
      end      
      if ~isempty(cplus)
         cJac(:,gcnt) = (cplus - cCurrent)/CHG(gcnt); 
      end           
   end
    xCurrent(gcnt) = temp;
end % for 








end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = nlsq(funfcn,x,verbosity,options,defaultopt,CostFunction,JAC,XDATA,YDATA,caller,varargin)
%NLSQ Helper function that solves non-linear least squares problems.
%   NLSQ is the core code for solving problems of the form:
%   min sum {FUN(X).^2} where FUN and X may be vectors or matrices.   
%             x
%

%   Copyright 1990-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/12/15 19:29:06 $

% ------------Initialization----------------
XOUT = x(:);
% numberOfVariables must be the name of this variable
numberOfVariables = length(XOUT);
msg = [];
how = [];
OUTPUT = [];
iter = 0;
EXITFLAG = 1;  %assume convergence
currstepsize = 0;

% Handle the output
if isfield(options,'OutputFcn')
    outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
else
    outputfcn = defaultopt.OutputFcn;
end
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end
stop = false;

% Handle the output
if isfield(options,'PlotFcns')
    plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
else
    plotfcns = defaultopt.PlotFcns;
end
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

formatstrFirstIter = ' %5.0f       %5.0f   %13.6g';
formatstr = ' %5.0f       %5.0f   %13.6g %12.3g    %12.3g';

% If lsqcurvefit, then add XDATA to objective's input list.
% xargin{1} will be updated at each iteration
if ~isempty(XDATA)
    xargin = {x,XDATA};
else
    xargin = {x};
end

% options
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
tolX = optimget(options,'TolX',defaultopt,'fast');
lineSearchType = strcmp(optimget(options,'LineSearchType',defaultopt,'fast'),'cubicpoly');
% If caller is fsolve, LevenbergMarquardt is not a member of the defaultopt structure,
% but options.LevenbergMarquardt will exist and will have been set to either 'on' or 'off'
% in fsolve. Thus, the call to optimget below will work as expected. 
levMarq = strcmp(optimget(options,'LevenbergMarquardt',defaultopt,'fast'),'on');
tolFun = optimget(options,'TolFun',defaultopt,'fast');
DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
typicalx = optimget(options,'TypicalX',defaultopt,'fast') ;
if ischar(typicalx)
   if isequal(lower(typicalx),'ones(numberofvariables,1)')
      typicalx = ones(numberOfVariables,1);
   else
      error('optim:nlsq:InvalidTypicalX','Option ''TypicalX'' must be a numerical value if not the default.')
   end
end
maxFunEvals = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxIter = optimget(options,'MaxIter',defaultopt,'fast');
if ischar(maxFunEvals)
   if isequal(lower(maxFunEvals),'100*numberofvariables')
      maxFunEvals = 100*numberOfVariables;
   else
      error('optim:nlsq:MaxFunEvals','Option ''MaxFunEvals'' must be a numeric value if not the default.')
   end
end

nfun=length(CostFunction);
numFunEvals = 1;
numGradEvals = 0;
MATX=zeros(3,1);
MATL=[CostFunction'*CostFunction;0;0];
FIRSTF=CostFunction'*CostFunction;
PCNT = 0;
EstSum=0.5;
% system of equations or overdetermined
if nfun >= numberOfVariables
   GradFactor = 0;  
else % underdetermined: singularity in JAC'*JAC or GRAD*GRAD' 
   GradFactor = 1;
end
done = false;
NewF = CostFunction'*CostFunction;


while ~done  
   % Work Out Gradients
   if ~(gradflag) || DerivativeCheck
      JACFD = zeros(nfun, numberOfVariables);  % set to correct size
      JACFD = systemSolve.finitedifferences(XOUT,x,funfcn,[],[],[],CostFunction,[],XDATA,YDATA,...
                              DiffMinChange,DiffMaxChange,typicalx,[],'all',[],[],[],[],[],[],...
                              false,varargin{:});       
      numFunEvals=numFunEvals+numberOfVariables;
      % In the particular case when there is only one function in more
      % than one variable, finitedifferences() returns JACFD in a column
      % vector, whereas nlsq() expects a single-row Jacobian:
      if (nfun == 1) && (numberOfVariables > 1)
         JACFD = JACFD';
      end
      % Gradient check
      if DerivativeCheck && gradflag         
         if isa(funfcn{3},'inline') 
            % if using inlines, the gradient is in funfcn{4}
            graderr(JACFD, JAC, formula(funfcn{4})); %
         else 
            % otherwise fun/grad in funfcn{3}
            graderr(JACFD, JAC,  funfcn{3});
         end
         DerivativeCheck = 0;
      else
         JAC = JACFD;
      end
   else
      x(:) = XOUT;
   end
   GradF = 2*(CostFunction'*JAC)'; %2*GRAD*CostFunction;
   
   %---------------Initialization of Search Direction------------------
   JacTJac = JAC'*JAC;
   if iter == 0
       OLDX = XOUT;
       OLDF = CostFunction;
       systemSolve.displayHeaders(verbosity,levMarq);
       
       % Initialize the output function.
       if haveoutputfcn || haveplotfcn
           [xOutputfcn, optimValues, stop] = systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'init',iter,numFunEvals, ...
               CostFunction,NewF,[],[],GradF,[],[],varargin{:});
           if stop
               [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
               return;
           end
       end
       
       if verbosity > 1
           disp(sprintf(formatstrFirstIter,iter,numFunEvals,NewF)); 
       end
       
       % 0th iteration
       if haveoutputfcn || haveplotfcn
           [xOutputfcn, optimValues, stop] = systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'iter',iter,numFunEvals, ...
               CostFunction,NewF,[],[],GradF,[],[],varargin{:});
           if stop
               [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
               return;
           end
       end

      % Disable the warnings about conditioning for singular and
      % nearly singular matrices
      warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
      warningstate2 = warning('off', 'MATLAB:singularMatrix');
      
      if condest(JacTJac)>1e16 
         SD=-(JacTJac+(norm(JAC)+1)*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)';
         if levMarq 
            GradFactor=norm(JAC)+1; 
         end
         how='COND';
      else
         SD=-(JacTJac+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)';
      end
      FIRSTF=NewF;
      OLDJ = JAC;
      GDOLD=GradF'*SD;
      % currstepsize controls the initial starting step-size.
      % If currstepsize has been set externally then it will
      % be non-zero, otherwise set to 1. Right now it's not
      % possible to set it externally.
      if currstepsize == 0, 
         currstepsize=1; 
      end
      LMorSwitchFromGNtoLM = false;
      if levMarq
         newf=JAC*SD+CostFunction;
         GradFactor=newf'*newf;
         SD=-(JacTJac+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)'; 
         LMorSwitchFromGNtoLM = true;
      end

      % Restore the warning states to their original settings
      warning(warningstate1)
      warning(warningstate2)

      newf=JAC*SD+CostFunction;
      EstSum=newf'*newf;
      status=0;
      if lineSearchType==0; 
         PCNT=1; 
      end
     
  else % iter >= 1
      gdnew=GradF'*SD;
      if verbosity > 1 
          num=sprintf(formatstr,iter,numFunEvals,NewF,currstepsize,gdnew);
          if LMorSwitchFromGNtoLM
              num=[num,sprintf('   %12.6g ',GradFactor)]; 
          end
          if isinf(verbosity)
              disp([num,'       ',how])
          else
              disp(num)
          end
      end
      
      if haveoutputfcn || haveplotfcn
          [xOutputfcn, optimValues, stop] = systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'iter',iter,numFunEvals, ...
              CostFunction,NewF,currstepsize,gdnew,GradF,SD,GradFactor,varargin{:});
          if stop
              [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
              return;
          end
      end
      
      %-------------Direction Update------------------


      if gdnew>0 && NewF>FIRSTF
         % Case 1: New function is bigger than last and gradient w.r.t. SD -ve
         % ... interpolate. 
         how='inter';
         [stepsize]=cubici1(NewF,FIRSTF,gdnew,GDOLD,currstepsize);
         currstepsize=0.9*stepsize;
      elseif NewF<FIRSTF
         %  New function less than old fun. and OK for updating 
         %         .... update and calculate new direction. 
         [newstep,fbest] =systemSolve.cubici3(NewF,FIRSTF,gdnew,GDOLD,currstepsize);
         if fbest>NewF,
            fbest=0.9*NewF; 
         end 
         if gdnew<0
            how='incstep';
            if newstep<currstepsize,  
               newstep=(2*currstepsize+1e-4); how=[how,'IF']; 
            end
            currstepsize=abs(newstep);
         else 
            if currstepsize>0.9
               how='int_step';
               currstepsize=min([1,abs(newstep)]);
            end
         end
         % SET DIRECTION.      
         % Gauss-Newton Method    
         LMorSwitchFromGNtoLM = true ;

         % Disable the warnings about conditioning for singular and
         % nearly singular matrices
         warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
         warningstate2 = warning('off', 'MATLAB:singularMatrix');

         if ~levMarq
            if currstepsize>1e-8 && condest(JacTJac)<1e16
               SD=JAC\(JAC*XOUT-CostFunction)-XOUT;
               if SD'*GradF>eps,
                  how='ERROR- GN not descent direction';
               end
               LMorSwitchFromGNtoLM = false;
            else
               if verbosity > 1
                   if currstepsize <= 1e-8
                       fprintf('Step size too small - Switching to LM method.\n');
                   else % condest(JacTJac) >= 1e16
                       fprintf('Iteration matrix ill-conditioned - Switching to LM method.\n');
                   end
               end
               how='CHG2LM';
               levMarq = true;
               currstepsize=abs(currstepsize);               
            end
         end
         
         if (LMorSwitchFromGNtoLM)      
            % Levenberg_marquardt Method N.B. EstSum is the estimated sum of squares.
            %                                 GradFactor is the value of lambda.
            % Estimated Residual:
            if EstSum>fbest
               GradFactor=GradFactor/(1+currstepsize);
            else
               GradFactor=GradFactor+(fbest-EstSum)/(currstepsize+eps);
            end
            SD=-(JacTJac+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)'; 
            currstepsize=1; 
            estf=JAC*SD+CostFunction;
            EstSum=estf'*estf;
         end
         gdnew=GradF'*SD;

         % Restore the warning states to their original settings
         warning(warningstate1)
         warning(warningstate2)
         
         OLDX=XOUT;
         % Save Variables
         FIRSTF=NewF;
         OLDG=GradF;
         GDOLD=gdnew;    
         
         % If quadratic interpolation set PCNT
         if lineSearchType==0, 
            PCNT=1; MATX=zeros(3,1); MATL(1)=NewF; 
         end
      else 
         % Halve Step-length
         how='Red_Step';
         if NewF==FIRSTF,
            msg = sprintf('No improvement in search direction: Terminating.');
            done = true;
            EXITFLAG = -4;
         else
            currstepsize=currstepsize/8;
            if currstepsize<1e-8
               currstepsize=-currstepsize;
            end
         end
      end
   end % if iter==0 
   
   %----------End of Direction Update-------------------
   iter = iter + 1;

   if lineSearchType==0, 
      PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; 
   end
   % Check Termination 
   if (GradF'*SD) < tolFun && ...
           max(abs(GradF)) < 10*(tolFun+tolX)
       msg = sprintf(['Optimization terminated: directional derivative along\n' ... 
                          ' search direction less than TolFun and infinity-norm of\n' ...
                          ' gradient less than 10*(TolFun+TolX).']);
       done = true; EXITFLAG = 1;
   elseif max(abs(SD))< tolX 
       msg = sprintf('Optimization terminated: magnitude of search direction less than TolX.');     
       done = true; EXITFLAG = 4;

   elseif numFunEvals > maxFunEvals
      msg = sprintf(['Maximum number of function evaluations exceeded.',...
                 ' Increase OPTIONS.MaxFunEvals.']);
      done = true;
      EXITFLAG = 0;
   elseif iter > maxIter
      msg = sprintf(['Maximum number of iterations exceeded.', ...
                        ' Increase OPTIONS.MaxIter.']);
      done = true;
      EXITFLAG = 0;
   else % continue
      XOUT=OLDX+currstepsize*SD;
      % Line search using mixed polynomial interpolation and extrapolation.
      if PCNT~=0
         % initialize OX and OLDF 
         OX = XOUT; OLDF = CostFunction;
         while PCNT > 0 && numFunEvals <= maxFunEvals
            x(:) = XOUT; 
            xargin{1} = x; % update x in list of input arguments to objective
            CostFunction = feval(funfcn{3},xargin{:},varargin{:});
            if ~isempty(YDATA)
               CostFunction = CostFunction - YDATA;
            end
            CostFunction = CostFunction(:);
            numFunEvals=numFunEvals+1;
            NewF = CostFunction'*CostFunction;
            % <= used in case when no improvement found.
            if NewF <= OLDF'*OLDF, 
               OX = XOUT; 
               OLDF=CostFunction; 
            end
            [PCNT,MATL,MATX,steplen,NewF,how]=systemSolve.searchq(PCNT,NewF,OLDX,MATL,MATX,SD,GDOLD,currstepsize,how);
            currstepsize=steplen;
            XOUT=OLDX+steplen*SD;
            if NewF==FIRSTF,  
               PCNT=0; 
            end
            if haveoutputfcn % Call output functions via callAllOptimOutputFcns wrapper
                stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,'interrupt',varargin{:});
                if stop
                    [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
                    return;
                end
            end
         end % end while
         XOUT = OX;
         CostFunction=OLDF;
         if numFunEvals>maxFunEvals 
            msg = sprintf(['Maximum number of function evaluations exceeded.',...
                            ' Increase OPTIONS.MaxFunEvals.']);
            done = true; 
            EXITFLAG = 0;
         end
      end % if PCNT~=0
   
      % Evaluate objective
      x(:)=XOUT;
      xargin{1} = x; % update x in list of input arguments to objective
      switch funfcn{1}
          case 'fun'
              CostFunction = feval(funfcn{3},xargin{:},varargin{:});
              if ~isempty(YDATA)
                  CostFunction = CostFunction - YDATA;
              end
              CostFunction = CostFunction(:);
              nfun=length(CostFunction);
              % JAC will be updated when it is finite-differenced
          case 'fungrad'
              [CostFunction,JAC] = feval(funfcn{3},xargin{:},varargin{:});
              if ~isempty(YDATA)
                  CostFunction = CostFunction - YDATA;
              end
              CostFunction = CostFunction(:);
              numGradEvals=numGradEvals+1;
          case 'fun_then_grad'
              CostFunction = feval(funfcn{3},xargin{:},varargin{:});
              if ~isempty(YDATA)
                  CostFunction = CostFunction - YDATA;
              end
              CostFunction = CostFunction(:);
              JAC = feval(funfcn{4},xargin{:},varargin{:});
              numGradEvals=numGradEvals+1;
          otherwise
              error('optim:nlsq:InvalidCalltype','Undefined calltype in LSQNONLIN.')
      end
      numFunEvals=numFunEvals+1;
      NewF = CostFunction'*CostFunction;

   end % Convergence testing
end  % while

NewF = CostFunction'*CostFunction;
gdnew=GradF'*SD;

if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'iter',iter,numFunEvals, ...
        CostFunction,NewF,currstepsize,GDOLD,GradF,SD,GradFactor,varargin{:});
    if stop
        [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = systemSolve.cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
        return;
    end
end

x(:)=XOUT;

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize=currstepsize;
OUTPUT.cgiterations = [];
OUTPUT.firstorderopt = [];

if levMarq
   OUTPUT.algorithm='medium-scale: Levenberg-Marquardt, line-search';
else
   OUTPUT.algorithm='medium-scale: Gauss-Newton, line-search';
end

if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues] = systemSolve.callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'done',iter,numFunEvals, ...
        CostFunction,NewF,currstepsize,GDOLD,GradF,SD,GradFactor,varargin{:});
    % Optimization done, so ignore "stop"
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayHeaders(verbosity,levMarq)

if verbosity>1
   if ~levMarq
      if isinf(verbosity)
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative   Line-search']);
      else
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative ']);
      end
   else
      if isinf(verbosity)
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative   Lambda       Line-search']);
      else
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative    Lambda']);
      end
   end
   disp(header)
end
end
%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,x,xOutputfcn,state,iter,numFunEvals, ...
    CostFunction,NewF,currstepsize,gdnew,GradF,SD,GradFactor,varargin)
% systemSolve.callOutputAndPlotFcns assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter', or 'done'. 

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.iteration = iter;
optimValues.funccount = numFunEvals;
optimValues.stepsize = currstepsize;
optimValues.directionalderivative = gdnew;
optimValues.gradient = GradF;
optimValues.firstorderopt = norm(GradF,Inf);
optimValues.searchdirection = SD;
optimValues.lambda = GradFactor;
if isequal(caller,'fsolve') 
   optimValues.fval = CostFunction; 
else % lsqnonlin, lsqcurvefit 
   optimValues.residual = CostFunction; 
   optimValues.resnorm = NewF; 
end 
xOutputfcn(:) = x;  % Set x to have user expected size
stop = false;
% Call output function
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        case 'interrupt'
            % No 'interrupt' case in nlsq
        otherwise
            error('optim:nlsq:UnknownStateInsystemSolve.callOutputAndPlotFcns','Unknown state in systemSolve.callOutputAndPlotFcns.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('optim:nlsq:UnknownStateInsystemSolve.callOutputAndPlotFcns','Unknown state in systemSolve.callOutputAndPlotFcns.')
    end
end
end
%--------------------------------------------------------------------------
function [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller)

x = xOutputfcn;
% CostFunction can be either 'fval' (fsolve) or 'residual'
if isequal(caller,'fsolve') 
    CostFunction = optimValues.fval;
else 
    CostFunction = optimValues.residual;
end
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.stepsize = optimValues.stepsize;
OUTPUT.cgiterations = [];
OUTPUT.firstorderopt = [];

if levMarq
   OUTPUT.algorithm='medium-scale: Levenberg-Marquardt, line-search';
else
   OUTPUT.algorithm='medium-scale: Gauss-Newton, line-search';
end

JAC = []; % May be in an inconsistent state
msg = 'Optimization terminated prematurely by user.';
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function step = quadi(x,f)
%QUADI Determine optimizer step from three points.
%   STEP = QUADI(x,f)
%   Finds the quadratic p(x) with p(x(1:3)) = f(1:3).
%   Returns the minimizer (or maximizer) of p(x).

%   Copyright 1990-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:38 $

% p(x) = a*x^2 + b*x + c.
% Solve [x.^2 x ones(3,1)]*[a b c]' = f.
% Minimum at p'(s) = 0,
% s = -b/(2*a) = (x1^2*(f2-f3) + x2^2*(f3-f1) + x3^2*(f1-f2))/ ...
%                 (x1*(f2-f3) + x2*(f3-f1) + x3*(f1-f2))/2

x = x(:); 
f = f(:);
g = f([2 3 1]) - f([3 1 2]);
step = ((x.*x)'*g)/(x'*g)/2;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pcnt, matl,matx,stepsize,fnew,how]=searchq(pcnt,fnew,oldx,matl,matx,sd,gdold,stepsize,how)
%SEARCHQ Line search routine for LSQNONLIN, and LSQCURVEFIT functions.
%   Performs line search procedure for least squares optimization. 
%   Uses Quadratic Interpolation. When finished pcnt returns 0.

%   Copyright 1990-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:40 $

if pcnt==1
% Case 1: Next point less than initial point. 
%     Increase step-length based on last gradient evaluation
    if fnew<matl(1)
% Quadratic Extrapolation using gradient of first point and 
% values of two other points.
        matl(2)=fnew;
        matx(2)=stepsize;
        newstep=-0.5*gdold*stepsize*stepsize/(fnew-gdold*stepsize-matl(1)+eps);
        if newstep<stepsize,how=[how,'QEF ']; newstep=1.2*stepsize; end
        stepsize=1.2*newstep;
        pcnt=2;
    else
% Case 2: New point greater than initial point. Decrease step-length.
        matl(3)=fnew;
        matx(3)=stepsize;
%Interpolate to get stepsize
        stepsize=max([1e-8*stepsize,-gdold*0.5*stepsize^2/(fnew-gdold*stepsize-matl(1)+eps)]);
        how=[how,'r'];
        pcnt=3;
    end
% Case 3: Last run was Case 1 (pcnt=2) and new point less than 
%     both of other 2. Replace. 
elseif pcnt==2  && fnew< matl(2)
    newstep=systemSolve.cubici2(gdold,[matl(1);matl(2);fnew],[matx(1);matx(2);stepsize]);
    if newstep<stepsize,how=[how, 'CEF ']; end
        matl(1)=matl(2);
        matx(1)=matx(2);
        matl(2)=fnew;
        matx(2)=stepsize;
        stepsize=min([newstep,1])+1.5*stepsize;
        stepsize=max([1.2*newstep,1.2*stepsize]);
        how=[how,'i'];
% Case 4: Last run was Case 2: (pcnt=3) and new function still 
%     greater than initial value.
elseif pcnt==3 && fnew>=matl(1)
    matl(2)=fnew;
    matx(2)=stepsize;
    if stepsize<1e-6
        newstep=-stepsize/2;
        % Give up if the step-size gets too small 
        % Stops infinite loops if no improvement is possible.
        if abs(newstep) < (eps * eps), pcnt = 0; end
    else
        newstep=systemSolve.cubici2(gdold,[matl(1);matl(3);fnew],[matx(1);matx(3);stepsize]);
    end
    matx(3)=stepsize;
    if isnan(newstep), stepsize=stepsize/2; else stepsize=newstep; end
    matl(3)=fnew;
    how=[how,'R'];
% Otherwise must have Bracketed Minimum so do quadratic interpolation.
%  ... having just increased step.
elseif pcnt==2 && fnew>=matl(2)
    matx(3)=stepsize;
    matl(3)=fnew;
    [stepsize]=systemSolve.cubici2(gdold,matl,matx);
    pcnt=4;
% ...  having just reduced step.
elseif  pcnt==3  && fnew<matl(1)
    matx(2)=stepsize;
    matl(2)=fnew;
    [stepsize]=systemSolve.cubici2(gdold,matl,matx);
    pcnt=4;
% Have just interpolated - Check to see whether function is any better 
% - if not replace.
elseif pcnt==4 
    pcnt=0;
    stepsize=abs(stepsize);
% If interpolation failed use old point.
    if fnew>matl(2),
        fnew=matl(2);
        how='f';
        stepsize=matx(2);       
    end
end %if pcnt==1
end
 

    end
    
end

