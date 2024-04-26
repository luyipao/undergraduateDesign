clear;
addpath('.\..\src\functions');
% testDopingFunction
% testLegendreBaseFunction
% testL2Projection
% testPiecewiseL2Projection
% testGetLegendreBasisInfo
rmpath('.\..\src\functions');

f(1)
function y = f(ii)
global ii
g(1);
end

function y = g(x)
global ii
y = h(ii);
end


function y = h(x)
y = x;
end
