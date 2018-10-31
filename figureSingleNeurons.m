function figureSingleNeurons(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end


figure('Position', [100 100 1500 800])

colAct = [27,158,119]/256;
colInact = [217,95,2]/256;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVE EXEMPLARY SESSION

activeExample = 7;
beta = getBetas(activeExample);
ax5 = subplot(245);
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [-.025 0 0 0])
text(-0.2, 1.2, 'D', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
scatterBetas(beta, colAct)
title('Active session')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVE SESSIONS TOGETHER

datasetInd = find(l.coefVar < .6);
[beta, r2] = getBetas(datasetInd);
subplot(246)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
text(-0.2, 1.2, 'E', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
scatterBetas(beta, colAct)
title('Active sessions')
text(-0.06, -0.045, ['R^2 = ' num2str(nanmean(r2),2) '\pm' num2str(nanstd(r2),2)])
text(-0.06, -0.055, ['n = ' num2str(size(beta,1))])

subplot(247)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
if ~isempty(beta)
    polarhistogram(atan2(beta(:,2), beta(:,1)), 20, 'FaceColor',colAct,'FaceAlpha',.5, ...
        'EdgeColor',colAct)
end
text(-0.2, 1.2, 'F', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INACTIVE EXEMPLARY SESSION

inactiveExample = 1;
beta = getBetas(inactiveExample);
subplot(241)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [-.025 0 0 0])
text(-0.2, 1.2, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
scatterBetas(beta, colInact)
title('Inactive session')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INACTIVE SESSIONS TOGETHER

datasetInd = find(l.coefVar > 1.2);
[beta, r2] = getBetas(datasetInd);
subplot(242)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
text(-0.2, 1.2, 'B', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
scatterBetas(beta, colInact)
title('Inactive sessions')
text(-0.06, -0.045, ['R^2 = ' num2str(nanmean(r2),2) '\pm' num2str(nanstd(r2),2)])
text(-0.06, -0.055, ['n = ' num2str(size(beta,1))])

subplot(243)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
polarhistogram(atan2(beta(:,2), beta(:,1)), 20, 'FaceColor', colInact,'FaceAlpha',.5, ...
    'EdgeColor', colInact)
text(-0.2, 1.2, 'C', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACROSS DATASETS

for i = 1:length(l.datasets)
    beta = getBetas(i);
    if ~isempty(beta)
        imbalanceILD(i) = sum(beta(:,1)<0)/sum(~isnan(beta(:,1)));
        imbalanceABL(i) = sum(beta(:,2)>0)/sum(~isnan(beta(:,1)));
        posGain(i) = nanmean(beta(:,1).*beta(:,2)<0);
    else
        imbalanceILD(i) = nan;
        imbalanceABL(i) = nan;
        posGain(i) = nan;
    end
end

subplot(3,4,4)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
text(-0.5, 1.2, 'G', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
y = imbalanceILD;
myscatter(l,y)
xlabel('CV')
[r,p] = corr(l.coefVar(:), y(:));
text(0.2, 0.2, ['r=' num2str(r,2) ', p=' num2str(p,3)])
ylabel({'Fraction of neurons','with contralateral preference'})

subplot(3,4,8)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
text(-0.5, 1.2, 'H', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
y = imbalanceABL;
myscatter(l,y)
xlabel('CV')
[r,p] = corr(l.coefVar(:), y(:));
text(0.2, 0.2, ['r=' num2str(r,2) ', p=' num2str(p,3)])
ylabel({'Fraction of neurons','with loud preference'})

subplot(3,4,12)
set(gca, 'OuterPosition', get(gca,'OuterPosition') + [.025 0 0 0])
text(-0.5, 1.2, 'I', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
y = posGain;
myscatter(l,y)
xlabel('CV')
[r,p] = corr(l.coefVar(:), y(:));
text(0.2, 0.2, ['r=' num2str(r,2) ', p=' num2str(p,3)])
ylabel({'Fraction of neurons','with positive gain'})
set(gca, 'YTick', [0 0.5 1])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSETS WITH SINGLE NEURONS ON SUBPLOT 1

beta = getBetas(activeExample);
if ~isempty(beta)
    ind = [33 44 3 17 34 114 49 95];
    axes(ax5)
    yspan5 = ylim;
    xspan5 = xlim;
    pos5 = get(gca, 'Position');
    hold on
    scatter(beta(ind,1), beta(ind,2), [], [0 0 0]);
    col = parula12();
    
    figax = axes('Position', [0 0 1 1]);
    axis([0 1 0 1])
    hold on
    axis off
    
    for i=1:length(ind)
        if i<=4
            xpos = .03;
            ypos = .5-i*.1;
        else
            xpos = .27;
            ypos = .5-(i-4)*.1;
        end
        axes('Position', [xpos ypos .05 .05])
        X = squeeze(nanmean(sum(l.datasets{activeExample}(ind(i),:,:,l.time>0 & l.time<0.15,:),4) / 0.15, 5));
        baseline = squeeze(sum(l.datasets{activeExample}(ind(i),:,:,l.time>-0.15 & l.time<0,:),4) / 0.15);
        baseline = nanmean(baseline(:));
        xlim([-5 47])
        plot(xlim, baseline*[1 1], 'k')
        hold on
        for j=1:3
            x = (1:12) + (j-1)*14;
            for ild=1:12
                plot(x(ild), X(ild,j), 'o', 'MarkerSize', 5+j*2, 'MarkerFaceColor', col(ild,:), 'MarkerEdgeColor', 'w')
            end
        end
        xlim([-5 47])
        axis off
        yspan = ylim;
        
        axes(figax)
        if i<=4
            plot([xpos+.055 -.005 + pos5(1)+pos5(3)/2 + beta(ind(i),1)/diff(xspan5)*pos5(3)], [ypos+(baseline-yspan(1))/(yspan(2)-yspan(1))*.05 pos5(2)+pos5(4)/2 + beta(ind(i),2)/diff(yspan5)*pos5(4)*.85], 'k')
        else
            plot([xpos-.005 .005 + pos5(1)+pos5(3)/2 + beta(ind(i),1)/diff(xspan5)*pos5(3)], [ypos+(baseline-yspan(1))/(yspan(2)-yspan(1))*.05 pos5(2)+pos5(4)/2 + beta(ind(i),2)/diff(yspan5)*pos5(4)*.85], 'k')
        end
    end
end

axes(figax)
plot([.74 .74], [0.1 0.94], 'k')

% export
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureSingleNeurons.pdf','-dpdf','-r0')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [beta, r2] = getBetas(datasetInd)
        beta = [];
        r2 = [];
        
        ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
        abl = [20 40 60];
        abl = abl-40;

        for d = 1:length(datasetInd)
            X  = squeeze(sum(l.datasets{datasetInd(d)}(:,:,:,l.time>0 & l.time<0.15,:),4) / 0.15);
            
            for n = 1:size(X,1)
                fr = X(n,:,:,:); 
                ind = ~isnan(fr(:));
                ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
                ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
%                 b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)], 0.01);

                % UNCOMMENT THE LINE BELOW TO INCLUDE THE MULTIPLICATIVE TERM
                b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1) ildAll(ind).*ablAll(ind)]);

                beta = [beta; b(1:2)'];                
                
                y = squeeze(nanmean(X(n,:,:,:),4));
                y = (y - mean(fr(ind)))/std(fr(ind));
                x1 = [ild' ild' ild'];
                x2 = repmat(abl, [length(ild) 1]);
                y = y(:);
                x1 = x1(:);
                x2 = x2(:);
                % yhat = x1*b(1) + x2*b(2) + b(3);
                yhat = x1*b(1) + x2*b(2) + b(3) + b(4)*x1.*x2;
                r2 = [r2; 1-sum((y-yhat).^2)/sum((y-mean(y)).^2)];  
            end
        end
    end

    function scatterBetas(beta, col)
        axis([-1 1 -1 1]*0.065)
        axis square
        hold on
        plot([0 0],ylim, 'k')
        plot(xlim,[0 0], 'k')
        xlabel('ILD coefficient')
        ylabel('ABL coefficient')
        if nargin==1
            col = [0 0 0];
        end
        scatter(beta(:,1), beta(:,2), 10, col, 'MarkerFaceColor', col);
        set(gca, 'XTick', [-0.05 0 0.05], 'YTick', [-0.05 0 0.05])
    end
end
