function figureInteractionTerm(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

sessions = {find(l.coefVar>1.2), find(l.coefVar<.6)};
col = {[217,95,2]/256, [27,158,119]/256};
titles = {'Inactive sessions', 'Active sessions'};

figure('Position', [100 100 1200 700])

% INACT or ACT
for mode = 1:2       
    betas = [];
    
    % go over all sessions
    for k = 1:length(sessions{mode})
        X = l.datasets{sessions{mode}(k)}(:,:,:,l.time>-0.05 & l.time<0.15,:) / 0.01;
        
        for n = 1:size(X,1)
            fr = squeeze(sum(X(n,:,:,6:end,:), 4));
            fr = nanmean(fr, 3); % regression on the tuning curves
            ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
            abl = [20 40 60];
            abl = abl-40;
            ildAll = bsxfun(@times, ild', ones([1 3 size(fr,3)]));
            ablAll = bsxfun(@times, abl,  ones([12 1 size(fr,3)]));
            ind = ~isnan(fr(:));
            [b,~,~,~,stats] = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ildAll(ind).*ablAll(ind) ones(sum(ind), 1)]);
            if stats(3) < 0.01
                betas = [betas; b(1:3)'];
            end
        end
    end
    
    subplot(2,3,(mode-1)*3+1)
    axis([-1 1 -1 1]*0.1)
    axis square
    hold on
    plot([0 0],ylim, 'k')
    plot(xlim,[0 0], 'k')
    xlabel('$\beta_\mathrm{ILD}$','Interpreter','latex')
    ylabel('$\beta_\mathrm{ABL}$','Interpreter','latex')
    scatter(betas(:,1), betas(:,2), 2, col{mode}, 'MarkerFaceColor', col{mode});
    set(gca, 'XTick', [-.1 -0.05 0 0.05 .1], 'YTick', [-.1 -0.05 0 0.05 .1])
    title(titles{mode})
        
    subplot(2,3,(mode-1)*3+2)
    axis([-1 1 -.05 .05]*0.1)
    axis square
    hold on
    plot([0 0],ylim, 'k')
    plot(xlim,[0 0], 'k')
    xlabel('$\beta_\mathrm{ILD}$','Interpreter','latex')
    ylabel('$\beta_\mathrm{inter}$','Interpreter','latex')
    scatter(betas(:,1), betas(:,3), 2, col{mode}, 'MarkerFaceColor', col{mode});
    set(gca, 'XTick', [-.1 -0.05 0 0.05 .1], 'YTick', [-0.005 0 0.005], 'YTickLabels', [-0.005 0 0.005])
    title(titles{mode})
    if mode==2
        [r,p] = corr(betas(:,1), betas(:,3));
        display(['Correlation ' num2str(r,2) ', p=' num2str(p,3)])
    end
    
    subplot(2,3,(mode-1)*3+3)
    axis([-1 1 -.05 .05]*0.1)
    axis square
    hold on
    plot([0 0],ylim, 'k')
    plot(xlim,[0 0], 'k')
    xlabel('$\beta_\mathrm{ABL}$','Interpreter','latex')
    ylabel('$\beta_\mathrm{inter}$','Interpreter','latex')
    scatter(betas(:,2), betas(:,3), 2, col{mode}, 'MarkerFaceColor', col{mode});
    set(gca, 'XTick', [-.1 -0.05 0 0.05 .1], 'YTick', [-0.005 0 0.005], 'YTickLabels', [-0.005 0 0.005])
    title(titles{mode})
end

letters = 'ABCDEF';
for i = 1:6
    subplot(2,3,i)
    text(-0.2, 1.14, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

% export (twice on purpose)
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperUnits','Inches','PaperSize', pos(3:4))
print(h,'figures/figureInteractionTerm.pdf','-dpdf','-r0')
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperUnits','Inches','PaperSize', pos(3:4))
print(h,'figures/figureInteractionTerm.pdf','-dpdf','-r0')
