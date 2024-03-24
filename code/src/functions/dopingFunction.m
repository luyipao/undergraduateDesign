function y = dopingFunction(x)
% 定义两个辅助函数，用于平滑过度
sinFunc1 = @(x) 249e15 * sin(20*pi*(x-0.075)) + 251e15;
sinFunc2 = @(x) 249e15 * sin(20*pi*(x-0.475)) + 251e15;

% 判断x的范围，计算各段的掺杂函数值
dopingDensity = @(x) (x>=0 & x<=0.1)*5e17 ...
    + (x>0.1 & x<0.15) .* sinFunc1(x) ...
    + (x>=0.15 & x<=0.45) * 2e15 ...
    + (x>0.45 & x<0.5) .* sinFunc2(x) ...
    + (x>=0.5 & x<=0.6) * 5e17;

% 计算掺杂函数值
y = dopingDensity(x);
end
