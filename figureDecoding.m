function figureDecoding(l)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 50 1300 900])

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
letters = 'IJK';
for abl=1:3
    for k=1:2
        for sub = 1:3
            subplot(3, 4, 9+sub)%subplots(sub))
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
        subplot(3,4, 9+abl)%5 + abl + (k-1)*4)
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
        ylabel('R^2')
        
%         bb = regress(accuracy(:,abl,k), [l.coefVar(:) l.coefVar(:)*0+1]);            
%         plot([0 2.2], [0 1; 2.2 1] * bb, style{k}, 'Color', 'r')        
     end
end


% NEUROMETRIC CURVES FOR ACT AND INACT

myfun = @(w,x) (w(1)-w(2))./(1+exp(-w(3)*(x-w(4)))) + w(2);
w0 = [1 0 1 0];
jnd = @(w) ((w(4)-1/w(3)*log((w(1)-0.75)/(0.75-w(2)))) - ...
            (w(4)-1/w(3)*log((w(1)-0.25)/(0.25-w(2))))) / 2;

decoding_psych(:,1:6,:) = 1 - decoding_psych(:,1:6,:);

ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20]';

inds = {find(l.coefVar<0.6), find(l.coefVar>1.2)};
letters = 'BA';
for m = 1:2
    %subplot(3,4,(m-1)*4+1)
    axes('Position', [.05 .1 + (m-1)*.45 .23 .35])
    axis([-20 20 0 1])
    hold on
    plot(xlim, [.5 .5], 'Color', [.5 .5 .5])
    plot([0 0], ylim, 'Color', [.5 .5 .5])
    ylabel('Fraction of ipsi classifications')
    xlabel('ILD')
    if m==1
        title('Active sessions')
    else
        title('Inactive sessions')
    end
    text(-0.1, 1.15, letters(m), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
    
    for i=1:3
        D = squeeze(mean(decoding_psych(inds{m},:,:),1));
        scatter(ild, D(:,i), [], 1-[i/3 i/3 i/3], 'filled')
        
        w = nlinfit(ild,D(:,i),myfun,w0);
        x = -20:.1:20;
        y = myfun(w,x);
        plot(x,y, 'k', 'LineWidth', 1+i/2)
        
        plot([2,5], [0+0.15*i, 0+0.15*i], 'k', 'LineWidth', 1+i/2)
        text(6, 0+0.15*i, ['ABL = ' num2str(i*20) ' dB'])
        text(6, 0+0.15*i-0.05, ['JND = ' num2str(jnd(w),2) ' dB'])
    end
end


% INTEGRAL PERFOMANCE PER SESSION BASED ON NEUROMETRIC CURVES

decoding_psych_perabl(:,1:6,:) = 1 - decoding_psych_perabl(:,1:6,:);

decoding_psych_shuffle = mean(decoding_psych_shuffle,4);
decoding_psych_shuffle(:,1:6,:) = 1 - decoding_psych_shuffle(:,1:6,:);

decoding_psych_shuffleTraining = mean(decoding_psych_shuffleTraining,4);
decoding_psych_shuffleTraining(:,1:6,:) = 1 - decoding_psych_shuffleTraining(:,1:6,:);

decodings = cat(4, decoding_psych, decoding_psych_perabl, decoding_psych_shuffle, decoding_psych_shuffleTraining);

warning('error', 'stats:nlinfit:IllConditionedJacobian');
warning('error', 'stats:nlinfit:IterationLimitExceeded');
warning('error', 'stats:nlinfit:ModelConstantWRTParam');
warning('error', 'MATLAB:rankDeficientMatrix');
warnings = [0 0 0 0];
for analysisType = 1:4
    for d = 1:length(l.datasets)
        for a = 1:3
            D = decodings(d,:,a,analysisType);
            try
                w = nlinfit(ild, D(:), myfun, w0);
            catch
                w = [nan nan nan nan];
                warnings(analysisType) = warnings(analysisType) + 1;
            end
            x = -20:.1:20;
            y = myfun(w,x);
            y(x<0) = 1-y(x<0);
            integralPerformances(d,a,analysisType) = mean(y);
        end
    end
    
    abl = repmat([20 40 60], [length(l.coefVar) 1]);
    abl = abl(:);
    cv = repmat(l.coefVar(:), [1 3]);
    cv = cv(:);
    y = integralPerformances(:, :, analysisType);
    y = y(:);
    [integralBetas(:, analysisType), ci] = regress(y, [cv abl abl.*cv cv*0+1]);
    integralCIs(:, :, analysisType) = ci;
end

% display(integralBetas)
% display(integralCIs)
display(warnings)

display(['Inact means: ' num2str(mean(integralPerformances(inds{2},:,1)))])
mu = mean(mean(integralPerformances(inds{2},:,1)));
delta = integralPerformances(inds{2},:,2) - integralPerformances(inds{2},:,1);
p = signrank(delta(:));
display(['Separate ABL diff: ' num2str(nanmean(delta(:))) '+-' num2str(nanstd(delta(:))), ', p=', num2str(p)])
delta = integralPerformances(inds{2},:,3) - integralPerformances(inds{2},:,1);
p = signrank(delta(:));
display(['Shuffle diff     : ' num2str(nanmean(delta(:))) '+-' num2str(nanstd(delta(:))), ', p=', num2str(p)])
acc2inf = @(q) q*log(q) + (1-q)*log(q) + log(2);
di = (acc2inf(mu)-acc2inf(mu+nanmean(delta(:))))/acc2inf(mu);
display(['delta I / I = ' num2str(di)])
delta = integralPerformances(inds{2},:,4) - integralPerformances(inds{2},:,1);
p = signrank(delta(:));
display(['Training sh diff : ' num2str(nanmean(delta(:))) '+-' num2str(nanstd(delta(:))), ', p=', num2str(p)])
di = (acc2inf(mu)-acc2inf(mu+nanmean(delta(:))))/acc2inf(mu);
display(['delta I / I = ' num2str(di)])

% Plot
letters = 'CDE';
for abl = 1:3
    for pl = 1:3
        subplot(3,4,1+pl)
        hold on
        w = integralBetas(:, 1);
        col = [.7 .7 .7];
        if pl==abl
            col = [0 0 0];
        end
        plot([0 2.2], [0 abl*20 0 1; 2.2 abl*20 2.2*abl*20 1] * w, ...
            'Color', col, 'LineWidth', 1+abl/2);
        text(-0.2, 1.2, letters(pl), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
    end
end

cols = {[0 0 0], [102,166,30]/256, [230,171,2]/256, [231,41,138]/256};
letters = 'FGH';
for analysisType = 2:4
    subplot(3, 4, 4 + analysisType)
    hold on
    axis([0 2.2 .5 1])
    for abl = 1:3
        w = integralBetas(:, 1);
        plot([0 2.2], [0 abl*20 0 1; 2.2 abl*20 2.2*abl*20 1] * w, ...
            'Color', [.7 .7 .7], 'LineWidth', 1+abl/2);
        w = integralBetas(:, analysisType);
        plot([0 2.2], [0 abl*20 0 1; 2.2 abl*20 2.2*abl*20 1] * w, ...
            'Color', cols{analysisType}, 'LineWidth', 1+abl/2);
    end
    ylabel('Accuracy')
    xlabel('CV')
    text(-0.2, 1.2, letters(analysisType-1), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end
subplot(346)
title('ABL separately')
subplot(347)
title('Shuffle')
subplot(348)
title('Training shuffle')

for abl=1:3
    subplot(3,4,1+abl)
    myscatter(l, integralPerformances(:,abl), 'noexamples')
    hold on
    ylim([0.5 1])
    title(['Integral acc. ABL=' num2str(abl*20)])
    xlabel('CV')
    ylabel('Accuracy')
end    

% Model comparison for integral performance
abl = repmat([20 40 60], [length(l.coefVar) 1]);
abl = abl(:);
cv = repmat(l.coefVar(:), [1 3]);
cv = cv(:);
y = integralPerformances(:, :, 1);
y = y(:);
beta = regress(y, [cv abl abl.*cv cv*0+1]);
rss = sum(([cv abl abl.*cv cv*0+1]*beta-y).^2);
n = length(l.coefVar)*3;
bic = n*log(rss/n) + 4*log(n);
aic = n*log(rss/n) + 2*4;
rss_saturated = 0;
for a=1:3
    ww = regress(integralPerformances(:,a,1), [l.coefVar(:) l.coefVar(:)*0+1]);
    rss_saturated = rss_saturated + sum(([l.coefVar(:) l.coefVar(:)*0+1] * ww - integralPerformances(:,a,1)).^2);
end
bic_saturated = n*log(rss_saturated/n) + 6*log(n);
aic_saturated = n*log(rss_saturated/n) + 2*6;
display(['Integral performance. My model:  AIC = ' num2str(aic) ', BIC = ' num2str(bic)])
display(['Integral performance. Saturated: AIC = ' num2str(aic_saturated) ', BIC = ' num2str(bic_saturated)])



h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureDecoding.pdf','-dpdf','-r0')
