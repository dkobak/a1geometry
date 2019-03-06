function summaryScatters(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

fig = figure('Position', [100 100 1500 900]);

for k = 1:length(l.datasets)
    pause(0.01)
    
    subplot(4,6,k)
    title(l.fileName{k}(6:end-12),'interpreter','none')
    
    beta = getBetas(k, true, false);
    scatterBetas(beta, [0,0,0], true)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [beta, r2] = getBetas(datasetInd, ifzscore, iffilter)
        beta = [];
        r2 = [];
        
        ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
        abl = [20 40 60];
        abl = abl - 40;

        for d = 1:length(datasetInd)
            X = squeeze(sum(l.datasets{datasetInd(d)}(:,:,:,l.time>0 & l.time<0.05,:),4) / 0.15);
            X = nanmean(X,4); % comment this line out to perform the regression on single trials
            
            for n = 1:size(X,1)
                fr = X(n,:,:,:); 
                
                ind = ~isnan(fr(:));
                ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
                ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
                if ifzscore
                    [b,~,~,~,stats] = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1) ildAll(ind).*ablAll(ind)]);
%                     [b,~,~,~,stats] = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)]);
                else
                    [b,~,~,~,stats] = regress(fr(ind), [ildAll(ind) ablAll(ind) ones(sum(ind), 1) ildAll(ind).*ablAll(ind)]);
                end
                
                % just in case a neuron has all values 0
                if sum(fr(ind))==0
                    b = b * 0;
                    stats(1) = 1;
                    stats(3) = 1;
                end
                                
                % filter out all neurons with model p>0.01
                if iffilter && (stats(3) > 0.01)
                    b = b * nan;
                    stats(1) = stats(1) * nan;
                end
                
                beta = [beta; b(1:2)'];
                r2 = [r2; stats(1)];
                
%                 y = squeeze(nanmean(X(n,:,:,:),4));
%                 y = (y - mean(fr(ind)))/std(fr(ind));
%                 x1 = [ild' ild' ild'];
%                 x2 = repmat(abl, [length(ild) 1]);
%                 y = y(:);
%                 x1 = x1(:);
%                 x2 = x2(:);
%                 % yhat = x1*b(1) + x2*b(2) + b(3);
%                 yhat = x1*b(1) + x2*b(2) + b(3) + b(4)*x1.*x2;
%                 r2 = [r2; 1-sum((y-yhat).^2)/sum((y-mean(y)).^2)];  
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