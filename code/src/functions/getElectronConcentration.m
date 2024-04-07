%% input: electronConcentrationCoeff C; mesh; poly max degree n
% function type: [)
function [electronConcentration,priElectronConcentration,electronConcentrationCells] = getElectronConcentration(C,mesh,n)
N = length(mesh) - 1;
%% init values
func1Cells = cell(N,1);
localAverageElectronConcentration = zeros(N,1);
for j = 1:N
    func1 = @(x) 0*x;
    for i = 1:n+1
        [Pi,~] = legendreBaseFunction(i-1, mesh(j),mesh(j+1));
        func1 = @(x) func1(x) + C((j-1)*(n+1)+i) * Pi(x);
    end
    func1Cells{j} = func1;
    localAverageElectronConcentration(j) = gaussLegendre(func1Cells{j}, mesh(j),mesh(j+1)) / (mesh(j+1)-mesh(j));
end

%%minmode limiter
%% coeffs(64)
% h = max(mesh(2:end)-mesh(1:end-1));
% M = 2 / 3 * 249*(20*pi)^2 * 10e2 * h^2;
% electronConcentrationCoeffs = zeros(N*(n+1),1);
% for j = 1:N
%     func1 = func1Cells{j};
%     ur = func1(mesh(j+1)) - localAverageElectronConcentration(j);
%     ul = localAverageElectronConcentration(j) - func1(mesh(j));
%     urmod = localAverageElectronConcentration(j) + minmod([ul localAverageElectronConcentration(mod(j,N)+1)-localAverageElectronConcentration(j) localAverageElectronConcentration(j)-localAverageElectronConcentration(mod(j-2+N,N)+1)], M);
%     ulmod =  localAverageElectronConcentration(j) - minmod([ur localAverageElectronConcentration(mod(j,N)+1)-localAverageElectronConcentration(j) localAverageElectronConcentration(j)-localAverageElectronConcentration(mod(j-2+N,N)+1)], M);
%     electronConcentrationCoeffs((j-1)*(n+1) + 1) = (mesh(j+1)-mesh(j)) * localAverageElectronConcentration(j) / sqrt(mesh(j+1)-mesh(j));
%     electronConcentrationCoeffs((j-1)*(n+1) + 2) = (urmod - ulmod) / (2*sqrt( 3/(mesh(j+1)-mesh(j)) ));
%     electronConcentrationCoeffs((j-1)*(n+1) + 3) = urmod - electronConcentrationCoeffs((j-1)*(n+1) + 1) / sqrt(mesh(j+1)-mesh(j)) - electronConcentrationCoeffs((j-1)*(n+1) + 2)*sqrt( 3/(mesh(j+1)-mesh(j)) );
%     electronConcentrationCoeffs((j-1)*(n+1) + 3) = electronConcentrationCoeffs((j-1)*(n+1) + 3) / sqrt( 5/(mesh(j+1)-mesh(j)) );
%     % electronConcentrationCoeffs((j-1)*(n+1) + 4) = 0;
% end

%% output
electronConcentration = @(x) 0*x;
priElectronConcentration = @(x) 0*x;
electronConcentrationCells = cell(N,1);
electronConcentrationCoeffs = C;
flux = 0;
for j = 1:N
    syms x
    func1 = @(x) 0*x;
    for i = 1:n+1
        [Pi,~] = legendreBaseFunction(i-1, mesh(j),mesh(j+1));
        func1 = @(x) func1(x) + electronConcentrationCoeffs((j-1)*(n+1)+i) * Pi(x);
    end
    func2 = matlabFunction(int(sym(func1)),'vars', x);
    flux = flux - func2(mesh(j));
    if j==1
        flux = 0;
    end
    electronConcentrationCells{j} = func1;
    if j==N
        electronConcentration = @(x) electronConcentration(x) + (x>=mesh(j) & x <= mesh(j+1)) .* func1(x);
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=mesh(j) & x <= mesh(j+1)) .* (func2(x) + flux);
    else
        electronConcentration = @(x) electronConcentration(x) + (x>=mesh(j) & x < mesh(j+1)) .* func1(x);
        priElectronConcentration = @(x) priElectronConcentration(x) + (x>=mesh(j) & x < mesh(j+1)) .* (func2(x) + flux);
    end
    flux = flux + func2(mesh(j+1));
end
end

function y = minmod(v,M)
% if abs(v(1)) <= M
%     y = v(1);
% else
%     y = all(sign(v) == sign(v(1))) * sign(v(1)) * min(abs(v));
% end
 y = all(sign(v) == sign(v(1))) * sign(v(1)) * min(abs(v));
end