function summaryPSTH(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

fig = figure('Position', [100 100 1500 900]);

for k = 1:length(l.datasets)
    pause(0.01)
    
    subplot(4,6,k)
    title(l.fileName{k}(6:end-12),'interpreter','none')
    
    X = l.datasets{k}(:,:,:,l.time>-0.05 & l.time<0.15,:) / 0.01;
    X = mean(X,1);
    X = nanmean(X,5);
    mua = squeeze(X);
    
    baseline = mua(:,:,1:5);
    baseline = mean(baseline(:));
    
    hold on
    axis([-0.05 0.15 0 max(mua(:))])
    plot(xlim, baseline*[1 1], 'Color', [.5 .5 .5])
    plot([0 0], ylim, 'Color', [.5 .5 .5])
    t = l.time(l.time>-0.05 & l.time<0.15);
    
    col = parula12();
    for ii=1:12
        for jj=1:3
            % y = squeeze(mua(ii,jj,:));
            % plot(t, y, 'Color', col(ii,:), 'LineWidth', jj);
            
            ntrial = find(~isnan(l.datasets{k}(1,ii,jj,1,:)), 1, 'last');
            yy = l.spikes{k}(:,ii,jj,:,1:ntrial);
            yy = yy(:);
            yy = yy(~isnan(yy));
            yy = histcounts(yy, [-0.3:0.001:0.3]);
            yy = yy / size(l.spikes{k},1) / ntrial;
            gaussKernel = 1/sqrt(2*pi)/5 * exp(-[-30:30].^2/5^2/2);
            yy = conv(yy, gaussKernel, 'same') * 1000;
            yy = yy(251:450);
            tt = (-0.05:0.001:0.15-.001) + 0.001/2;
            plot(tt,yy, 'Color', col(ii,:), 'LineWidth', jj);
        end
    end

end
