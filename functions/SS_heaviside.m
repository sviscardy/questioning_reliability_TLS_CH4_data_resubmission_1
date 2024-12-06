% y = SS_heaviside(x)
%
% Heaviside step function:
%
%   y = 0   for x<0
%   y = 1   for x>=0
%
% Author: sebastien.viscardy@aeronomie.be (March 29, 2023)
%%
function y = SS_heaviside(x)

nx      = length(x);
y       = zeros(1,nx);
y(x>=0) = 1;

end
