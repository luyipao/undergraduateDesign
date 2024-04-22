%
% addpath('C:\Users\12192\OneDrive - zju.edu.cn\documents\homeworks\undergraduateDesign\code\src\functions');
% setParameters
% a = 0;
% b = 0.6;
% n = 2;
% N = 100;
%
% mesh = Mesh(a,b,N,@(x) dopingFunction(x),n);
% mesh = mesh.DDModelDGFunction();

% f = @(x) sin(x);
% mesh = [0 pi/2 pi 2*pi];
% gaussLegendre(f,mesh(1:end-1),mesh(2:end))
% gaussLegendre(f,0,pi)
a = {1, 2, 3, 4};
b = {'A', 'B', 'C', 'D'};
temp =[ 1 0 1 0 ]; % '1'表示为真，'0'表示为假

c = a; % 创建一个和a一样的c
idx = find(temp); % 找到temp中为1的元素的索引位置
c = cellfun(@(x,y,index) y{i} 
c{idx} = b{idx}; % 使用comma-separated list assignment将b中对应位置的元素赋值给c。 

% 输出c的内容以验证结果
disp(c);