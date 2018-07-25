function figureMultiplicative(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

sessions = {find(l.coefVar>1.2), find(l.coefVar<.6)};

figure('Position', [100 100 1200 225])

% INACT or ACT
for mode = 1:2       
    betas = [];
    
    % go over all sessions
    for k = 1:length(sessions{mode})
        X = l.datasets{sessions{mode}(k)}(:,:,:,l.time>-0.05 & l.time<0.15,:) / 0.01;
        
        for n = 1:size(X,1)
            fr = squeeze(sum(X(n,:,:,6:end,:), 4));
            ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
            abl = [20 40 60];
            abl = abl-40;
            ildAll = bsxfun(@times, ild', ones([1 3 size(fr,3)]));
            ablAll = bsxfun(@times, abl,  ones([12 1 size(fr,3)]));
            ind = ~isnan(fr(:));
            b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ildAll(ind).*ablAll(ind) ones(sum(ind), 1)]);
            betas = [betas; b(1:3)'];
        end
    end
    
    if mode == 1
        subplot(141)
        [n,edges] = histcounts(betas(:,1),50);
        stairs(edges(1:end-1), n, 'Color', [217,95,2]/256)
        box off
        title('Inactive sessions')
        xlabel('ILD coefficient')
        xlim([-.05 .05])
        ylim([0 160])
        subplot(142)
        [n,edges] = histcounts(betas(:,2),50);
        stairs(edges(1:end-1), n, 'Color', [217,95,2]/256)
        box off        
        title('Inactive sessions')
        xlabel('ABL coefficient')
        xlim([-.05 .05])
        ylim([0 160])
        subplot(143)
        [n,edges] = histcounts(betas(:,3),50);
        stairs(edges(1:end-1)*1000, n, 'Color', [217,95,2]/256)
        box off
        title('Inactive sessions')
        xlabel('ILD*ABL coefficient \times{}10^3')
        xlim([-2 2])
        ylim([0 160])
    end
    
    if mode == 2
        subplot(144)
        scatter(betas(:,1), betas(:,3)*1000, 10, [27,158,119]/256, 'MarkerFaceColor', [27,158,119]/256)
        title('Active sessions')
        xlabel('ILD coefficient')
        ylabel('ILD*ABL coefficient \times{}10^3')
        xlim([-.05 .05])
        ylim([-2 2])
        [rho,p] = corr(betas(:,1), betas(:,3));
    end
end

letters = 'ABCD';
for i = 1:4
    subplot(1,4,i)
    text(-0.2, 1.14, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

% export (twice on purpose)
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperUnits','Inches','PaperSize', pos(3:4))
print(h,'figures/figureMultiplicative.pdf','-dpdf','-r0')
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperUnits','Inches','PaperSize', pos(3:4))
print(h,'figures/figureMultiplicative.pdf','-dpdf','-r0')
