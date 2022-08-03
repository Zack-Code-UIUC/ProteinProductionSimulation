function r = exp2rnd(mu1, mu2, tcrit, varargin)
%EXP2RND Random arrays from exponential distribution with two rates. 
%   R = EXP2RND(MU1,MU2,TCRIT) returns an array of random numbers chosen from the
%   exponential distribution with mean parameter MU1 up to time TCRIT and mean 
%   parameter MU2 afterwards. MU1 or MU2 or both can be scalars, in
%   which case the larger of the two determines the size of R. If they
%   are arrays of equal size, R will share their size. 
%
%   R = EXP2RND(MU1,MU2,TCRIT,M,N,...) or R = EXPRND(MU1,MU2,TCRIT,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   See also EXPRND, EXPCDF, EXPFIT, EXPINV, EXPLIKE, EXPPDF, EXPSTAT, RANDOM.
%
%   Author: Bruno Beltran <brunobeltran0@gmail.com>

%   EXP2RND uses the inversion method.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.

%   Copyleft 2014. Bruno Beltran.


if nargin < 1
    error(message('stats:exprnd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(3, mu1, mu2, tcrit, varargin{:});
if err > 0
    error(message('stats:exprnd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
mu1(mu1 < 0) = NaN;
mu2(mu2 < 0) = NaN;

% Generate uniform random values, and apply the exponential inverse CDF.
% if X = uniform, we want
% -mu1*ln(1-X)                    if X < 1-exp(-mu1*tcrit)
% -mu2*ln(1-X) + tcrit*(1-mu2/m1) otherwise
r = rand(sizeOut);
before_idx = r < 1 - exp(-1./mu1 .* tcrit);
if isscalar(mu1)
    mu1before = mu1;
    mu1after = mu1;
else
    mu1before = mu1(before_idx);
    mu1after = mu1(~before_idx);
end
if isscalar(mu2)
    mu2after = mu2;
else
    mu2after = mu2(~before_idx);
end
if ~isscalar(tcrit)
    tcrit = tcrit(~before_idx);
end
r(before_idx) = -mu1before.*log(1-r(before_idx));
r(~before_idx) = -mu2after.*log(1-r(~before_idx)) + tcrit.*(1-mu2after./mu1after);

if any(r < 0)
    error('');
end

