function figureIntro(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 100 500 800])

% Stimuli

col = parula12();

ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
abl = [20 40 60];

subplot(3,2,5)
hold on
for i=1:12
    for j=1:3
        plot(ild(i)/2+abl(j), -ild(i)/2+abl(j), 'o', 'MarkerSize', 5+j*2, 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'w')
    end
end

axis([0 80 0 80])
axis square
xlabel('Ipsi level (dB SPL)')
ylabel('Contra level (dB SPL)')
text(-0.2, 1.2, 'C', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)

% Raster plots

files = {'data/AC02-EXP01_Block-3_dataSet.mat', 'data/AC02-EXP08_Block-1_dataSet.mat'};
times = [60 6100];

for session=1:2
    load(files{session})

    subplot(3,2,(session-1)*2+[1 2])
    hold on
    tint = times(session) + [0 20];

    count = [];
    for i = 1:length(spikes)
        count(i) = sum(spikes{i}>tint(1) & spikes{i}<tint(2));
    end
    [~, nn] = sort(count);

    for ni = 1:length(spikes)
        n = nn(ni);
    
        if ~isempty(find(spikes{n}>tint(1) & spikes{n}<tint(2), 1))
            plot(repmat(spikes{n}(spikes{n}>tint(1) & spikes{n}<tint(2))', [2 1]) , [ni ni+0.9]', 'k')
        end
    end
    ylim([1 length(spikes)+1])
    
    plot(tint(1)+[0 2], [2 2], 'k', 'LineWidth', 2)
    text(tint(1)+.5, 10, '2 s')
    
    set(gca, 'visible', 'off')
end

subplot(3,2,[1 2])
text(-0.08, 1.1, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(3,2,[3 4])
text(-0.08, 1.1, 'B', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)

% CV vs fraction of down states across sessions

subplot(3,2,6)
myscatter(l, l.inactLevel)
ylim([0 80])
xlabel('Coefficient of variation')
ylabel('Fraction of down states (%)')

[r,p] = corr(l.coefVar(:), l.inactLevel(:));
text(.1, 75, ['$$r=' num2str(r,2) '$$'], 'Interpreter', 'latex')
text(-0.2, 1.2, 'D', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)

% export

h = gcf();
h.Renderer='Painters';
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize', [pos(3), pos(4)])
print(h, 'figures/figureIntro.pdf','-dpdf','-r0')

