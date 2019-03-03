function figureShanks(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

datasetInd = find(l.coefVar < .6);

figure('Position', [100 100 1200 400])
for d = 1:length(datasetInd)
    subplot(1, length(datasetInd), d)
    ii = datasetInd(d);
    X  = l.datasets{ii}(:,:,:,l.time>0 & l.time<0.15,:);
    betas = zeros(size(X,1), 2);
    
    % compute betas
    for n = 1:size(X,1)
        fr = squeeze(sum(X(n,:,:,:,:), 4)) / 0.15;
        fr = nanmean(fr, 3); % regression on the tuning curves
        ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
        abl = [20 40 60];
        abl = abl-40;
        ildAll = bsxfun(@times, ild', ones([1 3 size(fr,3)]));
        ablAll = bsxfun(@times, abl,  ones([12 1 size(fr,3)]));
        ind = ~isnan(fr(:));
        [b,~,~,~,stats] = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)]);
        betas(n,:) = b(1:2);
        if stats(3)>0.01
            betas(n,:) = [nan; nan];
        end
    end
    shanks = l.shanks{ii}(~isnan(betas(:,1)));
    betas = betas(~isnan(betas(:,1)),:);
    
    hold on
    for i=1:8
        g = scatter(betas(shanks==i,1), betas(shanks==i,2), '.');
        if sum(shanks==i) >= 5
            h = plotEllipse(betas(shanks==i,1:2), .9);
            h.set('Color', g.CData)
        end
    end
    axis([-1 1 -1 1]*0.1)
    axis square
    xlabel('$\beta_\mathrm{ILD}$','Interpreter','latex')
    ylabel('$\beta_\mathrm{ABL}$','Interpreter','latex')
end

% export
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureShanks.pdf','-dpdf','-r0')

