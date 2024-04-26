clear;
%f = @(x) sin(x);
f = @(x) dopingFunction(x);
a = 0.1; b = 0.15; n = 5;


 % draw 
 X = linspace(a,b,1000);
 fProj = L2Projection(f,n,a,b);
 plot(X,f(X),'-',X,fProj(X),'--');
legend('exact','L2 Projection')
 