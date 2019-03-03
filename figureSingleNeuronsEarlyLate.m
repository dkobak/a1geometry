function figureSingleNeuronsEarlyLate(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 100 1500 700])

colAct = [27,158,119]/256;
colInact = [217,95,2]/256;

datasetInd = find(l.coefVar > 1.2);
beta = getBetas(datasetInd, l.time>0 & l.time<0.05, true);
subplot(241)
scatterBetas(beta, colInact, true)
title('Inactive sessions. 0–50 ms')

beta = getBetas(datasetInd, l.time>0.1 & l.time<0.15, true);
subplot(242)
scatterBetas(beta, colInact, true)
title('Inactive sessions. 100–150 ms')

beta = getBetas(datasetInd, l.time>0 & l.time<0.05, false);
subplot(243)
scatterBetas(beta, colInact, false)
title({'Inactive sessions. 0–50 ms','Without z-scoring'})

beta = getBetas(datasetInd, l.time>0.1 & l.time<0.15, false);
subplot(244)
scatterBetas(beta, colInact, false)
title({'Inactive sessions. 100–150 ms','Without z-scoring'})

datasetInd = find(l.coefVar < .6);
beta = getBetas(datasetInd, l.time>0 & l.time<0.05, true);
subplot(245)
scatterBetas(beta, colAct, true)
% inspect = [34,35,60,91,95];
% scatter(beta(inspect,1), beta(inspect,2), 20, [1 0 0], 'MarkerFaceColor', 'r')
title('Active sessions. 0–50 ms')

beta = getBetas(datasetInd, l.time>0.1 & l.time<0.15, true);
subplot(246)
scatterBetas(beta, colAct, true)
title('Active sessions. 100–150 ms')

beta = getBetas(datasetInd, l.time>0 & l.time<0.05, false);
subplot(247)
scatterBetas(beta, colAct, false)
title({'Active sessions. 0–50 ms','Without z-scoring'})

beta = getBetas(datasetInd, l.time>0.1 & l.time<0.15, false);
subplot(248)
scatterBetas(beta, colAct, false)
title({'Active sessions. 100–150 ms','Without z-scoring'})

letters = 'ABCDEFGH';
for i = 1:8
    subplot(2,4,i)
    text(-0.2, 1.2, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

% export
h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureSingleNeuronsEarlyLate.pdf','-dpdf','-r0')


    function beta = getBetas(datasetInd, timeinterval, ifzscore)
        beta = [];
        
        ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
        abl = [20 40 60] - 40;

        for d = 1:length(datasetInd)
            X = squeeze(mean(l.datasets{datasetInd(d)}(:,:,:,timeinterval,:),4) / 0.01);
            X = nanmean(X,4);
            
            for n = 1:size(X,1)
                fr = X(n,:,:,:); 
                ind = ~isnan(fr(:));
                ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
                ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
                if ifzscore
                    [b,~,~,~,stats] = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ildAll(ind).*ablAll(ind) ones(sum(ind), 1)]);
                else
                    [b,~,~,~,stats] = regress(fr(ind), [ildAll(ind) ablAll(ind) ildAll(ind).*ablAll(ind) ones(sum(ind), 1)]);
                end
                
                % filter out all neurons with model p>0.01
                if stats(3) > 0.01
                    b = b * nan;
                    stats(1) = stats(1) * nan;
                end

                beta = [beta; b(1:2)'];
            end
        end   
    end

    function scatterBetas(beta, col, ifzscore)
        if ifzscore
            axis([-1 1 -1 1]*.1)
            set(gca, 'XTick', [-.1 -0.05 0 0.05 .1], 'YTick', [-.1 -0.05 0 0.05 .1])            
        else
            axis([-1 1 -1 1]*.5)
            set(gca, 'XTick', [-.5 -0.25 0 0.25 .5], 'YTick', [-.5 -0.25 0 0.25 .5])
        end
        axis square
        hold on
        plot([0 0],ylim, 'k')
        plot(xlim,[0 0], 'k')
        xlabel('$\beta_\mathrm{ILD}$','Interpreter','latex')
        ylabel('$\beta_\mathrm{ABL}$','Interpreter','latex')
        if nargin==1
            col = [0 0 0];
        end
        scatter(beta(:,1), beta(:,2), 2, col, 'MarkerFaceColor', col);
    end
end