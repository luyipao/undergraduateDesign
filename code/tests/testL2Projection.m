clear;
% 确定待处理函数，区间及初始数量
f = @(x) sin(x);
a = 0; b = 2*pi; n = 3;

% 初始化误差与模数
L2Error = [];

% 在区间逐步细化的情况下计算f的L2投影，并计算投影误差
for ii = 0:9
    X = linspace(a, b/2^ii, 1000);
    fProj = L2Projection(f, n, a, b/2^ii);
    L2Error = [L2Error norm(fProj(X)-f(X))];
end

%计算误差的对数比以观察收敛速率
disp(['---------' ' n = ' num2str(n) ' -----------------']);
disp(['accuracy : ' num2str(log2(L2Error(1:end-1) ./ L2Error(2:end)))]);
 
 % draw 
 X = linspace(a,2*b,1000);
 fProj = L2Projection(f,n,a,b);
 plot(X,f(X),'-',X,fProj(X),'--');
legend('exact','L2 Projection')
 