function figureDecodingIpsiContra(l)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 50 1000 350])

load decoding_results.mat

% LINEAR DECODING FOR IPSI/CONTRA

accuracy = cat(3, decoding_abl_contra_linear, decoding_abl_ipsi_linear);

% Do regressions
y = accuracy(:,:,:);
cv = repmat(l.coefVar(:), [1 3 2]);
abl = repmat([20 40 60], [length(l.coefVar) 1 2]);
ifipsi = repmat(shiftdim([0 1],-1), [length(l.coefVar) 3 1]);

y = y(:);
cv = cv(:);
abl = abl(:);
ifipsi = ifipsi(:);

X = [cv abl.*cv ifipsi.*cv 1+0*cv ifipsi abl]; % design matrix
[beta, betaconf] = regress(y, X, 0.001);
% display(beta)
% display(betaconf)
rss = sum((X*beta-y).^2);
n = length(y);
bic = n*log(rss/n) + size(X,2)*log(n);
aic = n*log(rss/n) + size(X,2)*2;
rss_saturated = 0;
for abl=1:3
    for k=1:2
        b = regress(accuracy(:,abl,k), [l.coefVar(:) l.coefVar(:)*0+1]);
        rss_saturated = rss_saturated + sum(([l.coefVar(:) l.coefVar(:)*0+1] * b - accuracy(:,abl,k)).^2);
    end
end
nparam = 12;
bic_saturated = n*log(rss_saturated/n) + nparam*log(n);
aic_saturated = n*log(rss_saturated/n) + nparam*2;
display(['Linear decoding. My model:  AIC = ' num2str(aic) ', BIC = ' num2str(bic)])
display(['Linear decoding. Saturated: AIC = ' num2str(aic_saturated) ', BIC = ' num2str(bic_saturated)])

% Plot
style = {'-', '--'};

%subplots = [6 7 8 10 11 12];
letters = 'ABC';
for abl=1:3
    for k=1:2
        for sub = 1:3
            subplot(1, 3, sub)%subplots(sub))
            hold on
            color = [.7 .7 .7];
            if sub == abl%(k-1)*3 + abl
                color = [0 0 0];
            end
            plot([0 2.2], [0 0 0 1 k-1 abl; 2.2 2.2*abl*20 2.2*(k-1) 1 k-1 abl] * beta, ...
                 style{k}, 'Color', color, 'LineWidth', 1 + abl/2);
            text(-0.2, 1.2, letters(sub), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
         end
    end
end

for abl=1:3
    for k=1:2
        subplot(1, 3, abl)%5 + abl + (k-1)*4)
        %y = accuracy(:,abl,k);
        myscatter(l, accuracy(:,abl,2), 'noexamples', 'open')
        myscatter(l, accuracy(:,abl,1), 'noexamples')
%        scatter(l.coefVar(ind), accuracy(ind,abl,2), [], [1 1 1], 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1)
%        scatter(l.coefVar(ind), accuracy(ind,abl,1), [], [0 0 0], 'filled', 'MarkerEdgeColor', 'none')

        ylim([-0.1 1])
        xlim([0 2.2])
        %if k==1
            title(['Linear, ABL=' num2str(abl*20)])
        %else
        %    title(['IPSI LINEAR, ABL=' num2str(abl*20)])
        %end
        xlabel('CV')
        if abl==1
            ylabel('R^2')
        end
%         bb = regress(accuracy(:,abl,k), [l.coefVar(:) l.coefVar(:)*0+1]);            
%         plot([0 2.2], [0 1; 2.2 1] * bb, style{k}, 'Color', 'r')        
     end
end


h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureDecodingIpsiContra.pdf','-dpdf','-r0')
