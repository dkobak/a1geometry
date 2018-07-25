function myscatter(l, y, noexample, ifopen)

ds = {find(l.coefVar<0.6), find(l.coefVar>1.2)};
ii = setdiff(1:length(l.datasets), [ds{1}; ds{2}]);

axis([0 2.2 0 ceil(max(y))])
axis square
hold on

if nargin <= 3 || isempty(ifopen)
    scatter(l.coefVar(ii),    y(ii),    50, [117,112,179]/256, 'filled', 'MarkerEdgeColor', 'w')
    scatter(l.coefVar(ds{1}), y(ds{1}), 50, [27,158,119]/256,  'filled', 'MarkerEdgeColor', 'w')
    scatter(l.coefVar(ds{2}), y(ds{2}), 50, [217,95,2]/256,    'filled', 'MarkerEdgeColor', 'w')
else
    scatter(l.coefVar(ii),    y(ii),    50, [117,112,179]/256,   'LineWidth', 1)
    scatter(l.coefVar(ds{1}), y(ds{1}), 50, [27,158,119]/256,    'LineWidth', 1  )
    scatter(l.coefVar(ds{2}), y(ds{2}), 50, [217,95,2]/256,      'LineWidth', 1)
end

if nargin<=2 || isempty(noexample)
    scatter(l.coefVar([1 7]), y([1 7]), 80, 'k')
end