function figureNeuronZoo(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 100 800 800])
col = parula12();

exemplaryDatasets = [1 7];
% offsets = [0 7];

for example = 1:2
    d = exemplaryDatasets(example);
    X = l.datasets{d}(:,:,:,l.time>0 & l.time<0.15,:);
    X = squeeze(sum(sum(nanmean(sum(X,4) / 0.15, 5),2),3));
    [~,ind] = sort(X,'descend');
    %subplots = bsxfun(@plus, [1:7]' + offsets(example), [0:14:84]);
    for n = 1:50
        subplots(n + (example-1)*50) = subplot(10,10, n + (example-1)*50); %subplots(n))
        X = squeeze(nanmean(sum(l.datasets{d}(ind(n),:,:,l.time>0 & l.time<0.15,:),4) / 0.15, 5));
        baseline = squeeze(sum(l.datasets{d}(ind(n),:,:,l.time>-0.15 & l.time<0,:),4) / 0.15);
        baseline = nanmean(baseline(:));
        xlim([0 42])
        ylim([0 20])
        hold on
        plot(xlim, baseline*[1 1], 'k')
        plot([0 0], ylim, 'k')
        
        for j=1:3
            x = (1:12) + (j-1)*14;
            for i=1:12
                plot(x(i), X(i,j), 'o', 'MarkerSize', (5+j*2)/2, 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'None')%[.99 .99 .99], 'LineWidth', .05)
            end
        end
        axis off
    end
end

for n=1:100
    pos = get(subplots(n), 'Position');
    if n<=50
        pos(2) = pos(2) + .025;
    else
    	pos(2) = pos(2) - .025;
	end
    set(subplots(n), 'Position', pos);
end

axes('Position', [0, 0, 1, 1])
axis([0 1 0 1])
hold on
axis off
plot([0.1,0.93], [.52 .52], 'k', 'LineWidth', 2)

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureNeuronZoo.pdf','-dpdf','-r0', '-painters')
