function assume(~, ~)
%ASSUME Assume symbolic relationship.
%   ASSUME(COND) sets the condition COND to be true. When you
%   set an assumption, the toolbox replaces any previous
%   assumptions on the free variables in COND with the new
%   assumption.
%   Assumptions on symfuns are not allowed.
%
%   See also SYM/ASSUME.

%   Copyright 2014 The MathWorks, Inc.
    error(message('symbolic:sym:AssumptionOnFunctionsNotSupported'));           
end
