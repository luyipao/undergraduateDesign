
close all
x = linspace(0, 0.6);

% [f{1,1,1}, E{1,1,1}] = M{1}{1,1}.getFuncs;
% [f{2,2,1}, E{2,2,1}] = M{1}{2,2}.getFuncs;
% % cmp mesh num N
% fig1 = figure;
% plot(x, f{1,1,1}.solve(x), x, f{2,2,1}.solve(x), 'o', 'LineWidth', 1);
% lgd = legend('N=100', 'N=200', 'Location', 'north');
% setFormat(lgd, 'x', 'n');
% filename = sprintf('DDTVDRK3Degree2mu0.75Nn.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% fig1_position = get(fig1, 'Position');
% 
% 
% fig2 = figure('Position', fig1_position);
% plot(x, E{1,1,1}(x), x, E{2,2,1}(x), 'o', 'LineWidth', 1);
% lgd = legend('N=100', 'N=200', 'Location', 'north');
% setFormat(lgd, 'x', 'E');
% filename = sprintf('DDTVDRK3Degree2mu0.75NE.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% [f{1,1,2}, E{1,1,2}] = M{2}{1,1}.getFuncs;
% % [f{2,2,2}, E{2,2,2}] = M{2}{2,2}.getFuncs;
% 
% % cmp mobility mu
% fig3 = figure('Position', fig1_position);
% plot( x, f{1,1,2}.solve(x), x, f{1,1,1}.solve(x),  'o', 'LineWidth', 1);
% lgd = legend('\mu = \mu(n_d)', '\mu = 0.75', 'Location', 'north');
% setFormat(lgd, 'x', 'n');
% filename = sprintf('DDTVDRK3Degree2N100mun.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% % 第四个图
% fig4 = figure('Position', fig1_position);
% plot( x, E{1,1,2}(x), x, E{1,1,1}(x), 'o', 'LineWidth', 1);
% lgd = legend('\mu = \mu(n_d)', '\mu = 0.75', 'Location', 'north');
% setFormat(lgd, 'x', 'E');
% filename = sprintf('DDTVDRK3Degree2N100muE.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% % cmp physical model DD and HF
% [f{1,1,3}, E{1,1,3}] = M{3}{1,1}.getFuncs;
% % [f{2,2,3}, E{2,2,3}] = M{3}{2,2}.getFuncs;
% fig5 = figure;
% plot(x, f{1,1,3}.solve(x), x, f{1,1,1}.solve(x), 'o', 'LineWidth', 1);
% lgd = legend('HF', 'DD', 'Location', 'north');
% setFormat(lgd, 'x', 'n');
% filename = sprintf('DDTVDRK3Degree2N100mu0.75modeln.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% fig6 = figure;
% plot(x, E{1,1,3}(x), x, E{1,1,1}(x), 'o', 'LineWidth', 1);
% lgd = legend('HF', 'DD', 'Location', 'north');
% setFormat(lgd, 'x', 'E');
% filename = sprintf('DDTVDRK3Degree2N100mu0.75modelE.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% % slope limiter 
% [f{1,1,4}, E{1,1,4}] = M{4}{1,1}.getFuncs;
% fig7 = figure;
% plot(x, f{1,1,4}.solve(x), x, f{1,1,1}.solve(x), 'o', 'LineWidth', 1);
% lgd = legend('limiter', 'no limiter', 'Location', 'north');
% setFormat(lgd, 'x', 'n');
% filename = sprintf('DDTVDRK3Degree2N100mu0.75limitern.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');
% 
% fig8 = figure;
% plot(x, E{1,1,4}(x), x, E{1,1,1}(x), 'o', 'LineWidth', 1);
% lgd = legend('limiter', 'no limiter', 'Location', 'north');
% setFormat(lgd, 'x', 'E');
% filename = sprintf('DDTVDRK3Degree2N100mu0.75limiterE.pdf');
% filename = fullfile('..\..\论文\docs\figure',filename);
% exportgraphics(gcf, filename, 'ContentType', 'vector');

fig9 = figure;
scatter(1:length(TV{1,1,1}), TV{1,1,1},'.');
lgd = legend('N=100');
setFormat(lgd, 'nt', 'TV');
filename = sprintf('TVDRKN100.pdf');
filename = fullfile('..\..\论文\docs\figure',filename);
exportgraphics(gcf, filename, 'ContentType', 'vector');


fig10 = figure;
scatter(1:length(TV{2,2,1}), TV{2,2,1},'.')
lgd = legend('N=100');
setFormat(lgd, 'nt', 'TV');
filename = sprintf('TVDRKN200.pdf');
filename = fullfile('..\..\论文\docs\figure',filename);
exportgraphics(gcf, filename, 'ContentType', 'vector');
% clear all figures
close all
function setFormat(lgd, xLabelName, yLabelName)
lgd.FontSize = 11; % 
lgd.Location = 'north'; %
set(gca,'XMinorTick', 'on', 'YMinorTick', 'on', 'linewidth',1.5); % change axes width
set(gca,'yticklabel',get(gca,'ytick'));
set(gca,'xticklabel',get(gca,'xtick'));
xlabel(xLabelName);
ylabel(yLabelName);
end