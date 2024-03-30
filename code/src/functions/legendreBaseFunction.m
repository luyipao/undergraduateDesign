% Function: Legendre
% Input: Degree 'n' as a vector and interval [a,b]
% Input condition: b > a
% Output: function handle PN  DPN
% Output condition: Each cell contains the corresponding base function

function [PN, DPN] = legendreBaseFunction(n,a,b)
    P0 = @(x) 0.*x + 1;
    DP0 = @(x) 0.*x + 0;
    P1 = @(x) 1.*x;
    DP1 = @(x) 0.*x + 1;
    if n==0
        P = P0;
        DP = DP0;
    elseif n == 1
        P = P1;
        DP = DP1;
    else
        for ii = 2:n
            P = @(x) ((2*ii - 1) .* x .* P1(x) - (ii - 1) .* P0(x)) / ii;
            DP = @(x) ((2*ii - 1) .* P1(x) + (2*ii - 1) .* x .* DP1(x) - (ii - 1) .* DP0(x)) / ii;
            P0 = P1;
            DP0 = DP1;
            P1 = P;
            DP1 = DP;
        end
    end
    PN = @(x) sqrt((2*n + 1) / (b - a)) * P((2.*x - b - a) / (b - a));
    DPN = @(x) sqrt(4*(2*n + 1) / (b - a)^3) * DP((2.*x - b - a) / (b - a));
end
