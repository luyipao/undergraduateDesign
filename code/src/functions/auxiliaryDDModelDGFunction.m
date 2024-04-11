% input: mesh; poly max degree n; nFunc
% output: auxiliary fucntion q coeff;


function [auxq,auxiliaryVarCoeffVec] = auxiliaryDDModelDGFunction(mesh, nFunc, n)
%% initlize
 X = linspace(0,0.6,1000);
 h = X(2) - X(1);
global RELAXATION_PARAMETER THETA
N = length(mesh) - 1;
meshSize = max(abs(mesh(2:end)-mesh(1:end-1)));

%% AC=F, C is auxiliaryVarCoeffVec
% cal A
% 预先分配矩阵
A = zeros((n+1)*N,(n+1)*N);
for j = 1:N
    A_j = zeros(n+1, n+1);
    for k = 1:n+1
        [Pk, DPk] = legendreBaseFunction(k-1,mesh(j),mesh(j+1));
        for l = 1:n+1
            [Pl, DPl] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            A_j(k,l) = quadgk(@(x) Pk(x).*Pl(x), mesh(j), mesh(j+1));
        end
    end
    A((j-1)*(n+1)+1 : j*(n+1), (j-1)*(n+1)+1 : j*(n+1)) = A_j;
end
% cal F
F = zeros(N*(n+1),1);
for j = 1:N
    for l = 1:n+1
        [Pl,DPl] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
        index = (j-1)*(n+1) + l;
        F(index) = nFunc(mesh(j+1)) * Pl(mesh(j+1)) - nFunc(mesh(j)) * Pl(mesh(j)) - quadgk(@(x) nFunc(x) .* DPl(x), mesh(j), mesh(j+1));
        F(index) = F(index) * sqrt(RELAXATION_PARAMETER * THETA);
    end
end
% cal C
auxiliaryVarCoeffVec = A\F;

%% init values
% func1Cells = cell(N,1);
% localAverageAuxq = zeros(N,1);
% for j = 1:N
%     func1 = @(x) 0*x;
%     for i = 1:n+1
%         [Pi,~] = legendreBaseFunction(i-1, mesh(j),mesh(j+1));
%         func1 = @(x) func1(x) + auxiliaryVarCoeffVec((j-1)*(n+1)+i) * Pi(x);
%     end
%     func1Cells{j} = func1;
%     localAverageAuxq(j) = quadgk(func1Cells{j}, mesh(j),mesh(j+1)) / (mesh(j+1)-mesh(j));
% end
%%minmode limiter
% h = mesh(2) - mesh(1);
% M = 0.1392 * 2 / 3 * 249*20*pi*20*pi * 20 * pi * 10e2 * h^2;
% auxqCoeffs = zeros(N*(n+1),1);
% for j = 1:N
%     func1 = func1Cells{j};
%     ur = func1(mesh(j+1)) - localAverageAuxq(j);
%     ul = localAverageAuxq(j) - func1(mesh(j));
%     urmod = localAverageAuxq(j) + minmod([ur localAverageAuxq(mod(j,N)+1)-localAverageAuxq(j) localAverageAuxq(j)-localAverageAuxq(mod(j-2+N,N)+1)], M);
%     ulmod =  localAverageAuxq(j) - minmod([ul localAverageAuxq(mod(j,N)+1)-localAverageAuxq(j) localAverageAuxq(j)-localAverageAuxq(mod(j-2+N,N)+1)], M);
%     auxqCoeffs((j-1)*(n+1) + 1) = (mesh(j+1)-mesh(j)) * localAverageAuxq(j) / sqrt(mesh(j+1)-mesh(j));
%     auxqCoeffs((j-1)*(n+1) + 2) = (urmod - ulmod) / (2*sqrt( 3/(mesh(j+1)-mesh(j)) ));
%     auxqCoeffs((j-1)*(n+1) + 3) = urmod - auxqCoeffs((j-1)*(n+1) + 1) / sqrt(mesh(j+1)-mesh(j)) - auxqCoeffs((j-1)*(n+1) + 2)*sqrt( 3/(mesh(j+1)-mesh(j)) );
%     auxqCoeffs((j-1)*(n+1) + 3)  = auxqCoeffs((j-1)*(n+1) + 3) / sqrt( 5/(mesh(j+1)-mesh(j)) );
%     % electronConcentrationCoeffs((j-1)*(n+1) + 4) = 0;
% end
% auxiliaryVarCoeffVec = auxqCoeffs;

%% output auxiliary function q
auxq = @(x) 0 * x;

for j = 1:N
    if j == 1
        for l = 1:n+1
            [Pl,~] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            auxq = @(x) auxq(x) + auxiliaryVarCoeffVec((j-1)*(n+1)+l) * (x >= mesh(j) & x <= mesh(j+1)) .* Pl(x);
        end
    else
        for l = 1:n+1
            [Pl,~] = legendreBaseFunction(l-1,mesh(j),mesh(j+1));
            auxq = @(x) auxq(x) + auxiliaryVarCoeffVec((j-1)*(n+1)+l) *(x >= mesh(j) & x < mesh(j+1)) .* Pl(x);
        end
    end
end
end


function y = minmod(v,M)
if abs(v(1)) <= M
    y = v(1);
else
    y = all(sign(v) == sign(v(1))) * sign(v(1)) * min(abs(v));
end
end

