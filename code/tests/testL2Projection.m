clear;
% 确定待处理函数，区间及初始数量
%f = @(x) sin(x);
f = @(x) dopingFunction(x) / 10e17;
a = 0.096; b = 0.1020; n = 4;


 % draw 
 X = linspace(a,b,1000);
 fProj = L2Projection(f,n,a,b);
 plot(X,f(X),'-',X,fProj(X),'--');
legend('exact','L2 Projection')
 