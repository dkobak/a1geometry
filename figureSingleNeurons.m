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
subplot(245)
scatterBetas(beta, colAct)
title('Active session')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVE SESSIONS TOGETHER

datasetInd = find(l.coefVar < .6);
[beta, r2] = getBetas(datasetInd);
subplot(246)
scatterBetas(beta, colAct)
title('Active sessions')
text(-0.06, -0.045, ['R2 = ' num2str(nanmean(r2),2) ' +- ' num2str(nanstd(r2),2)])
text(-0.06, -0.055, ['n = ' num2str(size(beta,1))])

subplot(247)
if ~isempty(beta)
    polarhistogram(atan2(beta(:,2), beta(:,1)), 20, 'FaceColor',colAct,'FaceAlpha',.5, ...
        'EdgeColor',colAct)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INACTIVE EXEMPLARY SESSION

inactiveExample = 1;
beta = getBetas(inactiveExample);
subplot(241)
scatterBetas(beta, colInact)
title('Inactive session')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INACTIVE SESSIONS TOGETHER

datasetInd = find(l.coefVar > 1.2);
[beta, r2] = getBetas(datasetInd);
subplot(242)
scatterBetas(beta, colInact)
title('Inactive sessions')
text(-0.06, -0.045, ['R2 = ' num2str(nanmean(r2),2) ' +- ' num2str(nanstd(r2),2)])
text(-0.06, -0.055, ['n = ' num2str(size(beta,1))])

subplot(243)
polarhistogram(atan2(beta(:,2), beta(:,1)), 20, 'FaceColor', colInact,'FaceAlpha',.5, ...
    'EdgeColor', colInact)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACROSS DATASETS

for i = 1:length(l.datasets)
    beta = getBetas(i);
    if ~isempty(beta)
        imbalanceILD(i) = sum(beta(:,1)<0)/sum(~isnan(beta(:,1)));
        imbalanceABL(i) = sum(beta(:,2)>0)/sum(~isnan(beta(:,1)));
    else
        imbalanceILD(i) = nan;
        imbalanceABL(i) = nan;
    end
end

subplot(244)
y = imbalanceILD;
myscatter(l,y)
xlabel('CV')
[r,p] = corr(l.coefVar(:), y(:));
text(0.2, 0.2, ['r=' num2str(r,2) ', p=' num2str(p,3)])
ylabel('Fraction of contra neurons')

subplot(248)
y = imbalanceABL;
myscatter(l,y)
xlabel('CV')
[r,p] = corr(l.coefVar(:), y(:));
text(0.2, 0.2, ['r=' num2str(r,2) ', p=' num2str(p,3)])
ylabel('Fraction of loud neurons')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSETS WITH SINGLE NEURONS ON SUBPLOT 1

letters = 'ABCDEFGH';
subplots = [1 2 3 5 6 7 4 8];
for i = 1:8
    subplot(2,4,subplots(i))
    text(-0.2, 1.2, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

beta = getBetas(activeExample);
if ~isempty(beta)
ind = [33 44 3 17 34 114 49 95];
subplot(245)
hold on
scatter(beta(ind,1), beta(ind,2), [], [0 0 0]);
col = parula12();
for i=1:length(ind)
    if i<=4
        xpos = .08;
        ypos = .5-i*.1;
    else
        xpos = .27;
        ypos = .5-(i-4)*.1;
    end
    axes('Position', [xpos ypos .05 .05])
    X = squeeze(nanmean(sum(l.datasets{activeExample}(ind(i),:,:,l.time>0 & l.time<0.15,:),4) / 0.15, 5));
    baseline = squeeze(sum(l.datasets{activeExample}(ind(i),:,:,l.time>-0.15 & l.time<0,:),4) / 0.15);
    baseline = nanmean(baseline(:));        
    xlim([0 42])
    plot(xlim, baseline*[1 1], 'k')
    hold on
    for j=1:3
        x = (1:12) + (j-1)*14;
        for i=1:12
            plot(x(i), X(i,j), 'o', 'MarkerSize', 5+j*2, 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'w')
        end
    end
    axis off    
end
end

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
                b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)], 0.01);

                % UNCOMMENT THE LINE BELOW TO INCLUDE THE MULTIPLICATIVE TERM
                % b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1) ildAll(ind).*ablAll(ind)]);

                beta = [beta; b(1:2)'];                
                
                y = squeeze(nanmean(X(n,:,:,:),4));
                y = (y - mean(fr(ind)))/std(fr(ind));
                x1 = [ild' ild' ild'];
                x2 = repmat(abl, [length(ild) 1]);
                y = y(:);
                x1 = x1(:);
                x2 = x2(:);
                yhat = x1*b(1) + x2*b(2) + b(3);
                % yhat = x1*b(1) + x2*b(2) + b(3) + b(4)*x1.*x2;
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
        xlabel('ILD')
        ylabel('ABL')
        if nargin==1
            col = [0 0 0];
        end
        scatter(beta(:,1), beta(:,2), 10, col, 'MarkerFaceColor', col);
        set(gca, 'XTick', [-0.05 0 0.05], 'YTick', [-0.05 0 0.05])
    end
end
