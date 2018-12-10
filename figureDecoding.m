function figureDecoding(l)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 50 800 950])

load decoding_results_new.mat

% NEUROMETRIC CURVES FOR ACT AND INACT

myfun = @(w,x) (w(1)-w(2))./(1+exp(-w(3)*(x-w(4)))) + w(2);
w0 = [1 0 1 0];
jnd = @(w) ((w(4)-1/w(3)*log((w(1)-0.75)/(0.75-w(2)))) - ...
            (w(4)-1/w(3)*log((w(1)-0.25)/(0.25-w(2))))) / 2;

decoding_psych(:,1:6,:) = 1 - decoding_psych(:,1:6,:);

ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20]';

inds = {find(l.coefVar>1.2), find(l.coefVar<0.6)};
letters = 'AB';
for m = 1:2
    %subplot(3,4,(m-1)*4+1)
    %axes('Position', [.05 .1 + (m-1)*.45 .23 .35])
%     subplot(7,6,[1 2 3 7 8 9 13 14 15] + (m-1)*3)
    axes('Position', [.1 + (m-1)*.45 .65 .4 .3])

    axis([-20 20 0 1])
    hold on
    plot(xlim, [.5 .5], 'Color', [.5 .5 .5])
    plot([0 0], ylim, 'Color', [.5 .5 .5])
    xlabel('ILD')
    if m==1
        title('Inactive sessions')
        ylabel('Fraction of ipsi classifications')
    else
        title('Active sessions')
    end
    text(-0.1, 1.1, letters(m), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
    
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
% mu = mean(mean(integralPerformances(inds{2},:,1)));
delta = integralPerformances(inds{2},:,2) - integralPerformances(inds{2},:,1);
p = signrank(delta(:));
display(['Separate ABL diff: ' num2str(nanmean(delta(:))) '+-' num2str(nanstd(delta(:))), ', p=', num2str(p)])
delta = integralPerformances(inds{2},:,3) - integralPerformances(inds{2},:,1);
p = signrank(delta(:));
display(['Shuffle diff     : ' num2str(nanmean(delta(:))) '+-' num2str(nanstd(delta(:))), ', p=', num2str(p)])
% acc2inf = @(q) q*log(q) + (1-q)*log(q) + log(2);
% di = (acc2inf(mu)-acc2inf(mu+nanmean(delta(:))))/acc2inf(mu);
% display(['delta I / I = ' num2str(di)])
delta = integralPerformances(inds{2},:,4) - integralPerformances(inds{2},:,1);
p = signrank(delta(:));
display(['Training sh diff : ' num2str(nanmean(delta(:))) '+-' num2str(nanstd(delta(:))), ', p=', num2str(p)])
% di = (acc2inf(mu)-acc2inf(mu+nanmean(delta(:))))/acc2inf(mu);
% display(['delta I / I = ' num2str(di)])

% Plot
letters = 'CDE';
for abl = 1:3
    for pl = 1:3
        %subplot(3,4,1+pl)
%         subplot(7,6,[19 20 25 26] + (pl-1)*2)
        if abl==1
            ax(pl) = axes('Position', [.1 + (pl-1)*.3 .35 .25 .22]);
        else
            axes(ax(pl));
        end

        hold on
        w = integralBetas(:, 1);
        col = [.7 .7 .7];
        if pl==abl
            col = [0 0 0];
        end
        plot([0 2.2], [0 abl*20 0 1; 2.2 abl*20 2.2*abl*20 1] * w, ...
            'Color', col, 'LineWidth', 1+abl/2);
        text(-0.2, 1.15, letters(pl), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
    end
end

for abl=1:3
%     subplot(3,4,1+abl)
    axes(ax(abl))
    myscatter(l, integralPerformances(:,abl), 'noexamples')
    hold on
    ylim([0.5 1])
    title(['Integral acc., ABL=' num2str(abl*20)])
    xlabel('CV')
    if abl==1
        ylabel('Accuracy')
    end
end  


cols = {[0 0 0], [102,166,30]/256, [230,171,2]/256, [231,41,138]/256};
titles = {'','ABL separately', 'Shuffle', 'Training shuffle'};
letters = 'FGH';
for analysisType = 2:4
    %subplot(3, 4, 4 + analysisType)
%     subplot(7,6,[31 32 37 38] + (analysisType-2)*2)
    axes('Position', [.1 + (analysisType-2)*.3 .05 .25 .22])
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
    if analysisType==2
        ylabel('Accuracy')
    end
    xlabel('CV')
    text(-0.2, 1.15, letters(analysisType-1), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
    title(titles{analysisType})
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
