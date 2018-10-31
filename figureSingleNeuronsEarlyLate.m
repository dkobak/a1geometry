function figureSingleNeuronsEarlyLate(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 100 1200 700])

colAct = [27,158,119]/256;
colInact = [217,95,2]/256;

datasetInd = find(l.coefVar > 1.2);
beta = getBetas(datasetInd, l.time>0 & l.time<0.15);
subplot(231)
scatterBetas(beta, colInact)
title('Inactive sessions. 0–150 ms')

beta = getBetas(datasetInd, l.time>0 & l.time<0.05);
subplot(232)
scatterBetas(beta, colInact)
title('Inactive sessions. 0–50 ms')

beta = getBetas(datasetInd, l.time>0.1 & l.time<0.15);
subplot(233)
scatterBetas(beta, colInact)
title('Inactive sessions. 100–150 ms')

datasetInd = find(l.coefVar < .6);
beta = getBetas(datasetInd, l.time>0 & l.time<0.15);
subplot(234)
scatterBetas(beta, colAct)
title('Active sessions. 0–150 ms')

% beta = getBetas(datasetInd, l.time>0.02 & l.time<0.04);
beta = getBetas(datasetInd, l.time>0 & l.time<0.05);
subplot(235)
scatterBetas(beta, colAct)
title('Active sessions. 0–50 ms')

beta = getBetas(datasetInd, l.time>0.1 & l.time<0.15);
subplot(236)
scatterBetas(beta, colAct)
title('Active sessions. 100–150 ms')

letters = 'ABCDEF';
for i = 1:6
    subplot(2,3,i)
    text(-0.2, 1.2, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

% export
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureSingleNeuronsEarlyLate.pdf','-dpdf','-r0')


    function beta = getBetas(datasetInd, timeinterval)
        beta = [];
        
        ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
        abl = [20 40 60] - 40;

%         XX = [];
        for d = 1:length(datasetInd)
            X = squeeze(mean(l.datasets{datasetInd(d)}(:,:,:,timeinterval,:),4) / 0.01);
%             XX = cat(1,XX,nanmean(X,4)); 
            for n = 1:size(X,1)
                fr = X(n,:,:,:); 
                ind = ~isnan(fr(:));
                ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
                ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
                b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ildAll(ind).*ablAll(ind) ones(sum(ind), 1)]);
                beta = [beta; b(1:2)'];                
            end
        end
        
%         fr = mean(XX,1);
%         ind = ~isnan(fr(:));
%         ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
%         ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
%         b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ildAll(ind).*ablAll(ind) ones(sum(ind), 1)]);
%         beta = b(1:2)';    
    end

    function scatterBetas(beta, col)
        axis([-1 1 -1 1]*0.05)
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