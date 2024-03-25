clear; clc;
addpath('./functions');
% 假设你的table名为T
Age = [38;43;38;40;11];
Smoker = logical([1;0;1;0;1]);
T = table(Age, Smoker);


table2latex(T);