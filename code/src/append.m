% 定义设置格式的函数

x = linspace(0, 0.6);
% 获取函数和能量值
[f{1,1,1}, E{1,1,1}] = M{1}{1,1}.getFuncs;
[f{2,2,1}, E{2,2,1}] = M{1}{2,2}.getFuncs;

fig1 = figure;
plot(x, f{1,1,1}.solve(x), x, f{2,2,1}.solve(x), 'o', 'LineWidth', 1);
lgd = legend('N=100', 'N=200', 'Location', 'north');
setFormat(lgd, 'x', 'n');
filename = sprintf('DDTVDRK3Degree2mu0.75Nn.pdf');
filename = fullfile('..\..\论文\docs\figure',filename);
exportgraphics(gcf, filename, 'ContentType', 'vector');
fig1_position = get(fig1, 'Position');


fig2 = figure('Position', fig1_position);
plot(x, E{1,1,1}(x), x, E{2,2,1}(x), 'o', 'LineWidth', 1);
lgd = legend('N=100', 'N=200', 'Location', 'north');
setFormat(lgd, 'x', 'E');
filename = sprintf('DDTVDRK3Degree2mu0.75NE.pdf');
filename = fullfile('..\..\论文\docs\figure',filename);
exportgraphics(gcf, filename, 'ContentType', 'vector');

[f{1,1,2}, E{1,1,2}] = M{2}{1,1}.getFuncs;
[f{2,2,2}, E{2,2,2}] = M{2}{2,2}.getFuncs;

% 第三个图
fig3 = figure('Position', fig1_position);
plot(x, f{1,1,1}.solve(x), x, f{1,1,2}.solve(x), 'o', 'LineWidth', 1);
lgd = legend('\mu = 0.75', '\mu = \mu(n_d)', 'Location', 'north');
setFormat(lgd, 'x', 'n');
filename = sprintf('DDTVDRK3Degree2N100mun.pdf');
filename = fullfile('..\..\论文\docs\figure',filename);
exportgraphics(gcf, filename, 'ContentType', 'vector');
% 第四个图
fig4 = figure('Position', fig1_position);
plot(x, E{1,1,1}(x), x, E{1,1,2}(x), 'o', 'LineWidth', 1);
lgd = legend('\mu = 0.75', '\mu = \mu(n_d)', 'Location', 'north');
setFormat(lgd, 'x', 'E');
filename = sprintf('DDTVDRK3Degree2N100muE.pdf');
filename = fullfile('..\..\论文\docs\figure',filename);
exportgraphics(gcf, filename, 'ContentType', 'vector');

function setFormat(lgd, xLabelName, yLabelName)
lgd.FontSize = 11; % 
lgd.Location = 'north'; %
set(gca,'XMinorTick', 'on', 'YMinorTick', 'on', 'linewidth',1.5); % change axes width
xlabel(xLabelName);
ylabel(yLabelName);
end