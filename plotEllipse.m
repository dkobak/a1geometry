function h = plotEllipse(data, coverage)

if nargin==1
    coverage = .5;
end

scale = sqrt(chi2inv(coverage, 2));
%scale = 1;

data = data(~isnan(sum(data,2)), 1:2);
data0 = bsxfun(@minus, data, mean(data));

[v, d] = eig(cov(data0));
[d, ind] = sort(diag(d), 'descend');
d = diag(d);
v = v(:, ind);

t = linspace(0,2*pi,100);
e = [cos(t) ; sin(t)];
vv = v*sqrt(d) * scale;               
e = bsxfun(@plus, vv*e, mean(data)');

h = plot(e(1,:), e(2,:));