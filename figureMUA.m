function figureMUA(l)

% This function produces MUA Figure

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 100 1500 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUA PSTH (time-resolved)

% Two example sessions to use (as numbered in l structure array)
ds = [1 7];

for k=1:2    
    X = l.datasets{ds(k)}(:,:,:,l.time>-0.05 & l.time<0.15,:) / 0.01;
    X = mean(X,1);
    X = nanmean(X,5);
    mua = squeeze(X);
    
    baseline = mua(:,:,1:5);
    baseline = mean(baseline(:));
    
    subplot(2,4,(k-1)*4+1)
    hold on
    axis([-0.05 0.15 0 5])
    plot(xlim, baseline*[1 1], 'Color', [.5 .5 .5])
    plot([0 0], ylim, 'Color', [.5 .5 .5])
    t = l.time(l.time>-0.05 & l.time<0.15);
    
    col = parula12();
    for ii=1:12
        for jj=1:3
            % y = squeeze(mua(ii,jj,:));
            % plot(t, y, 'Color', col(ii,:), 'LineWidth', jj);
            
            ntrial = find(~isnan(l.datasets{ds(k)}(1,ii,jj,1,:)), 1, 'last');
            yy = l.spikes{ds(k)}(:,ii,jj,:,1:ntrial);
            yy = yy(:);
            yy = yy(~isnan(yy));
            yy = histcounts(yy, [-0.3:0.001:0.3]);
            yy = yy / size(l.spikes{ds(k)},1) / ntrial;
            gaussKernel = 1/sqrt(2*pi)/5 * exp(-[-30:30].^2/5^2/2);
            yy = conv(yy, gaussKernel, 'same') * 1000;
            yy = yy(251:450);
            tt = (-0.05:0.001:0.15-.001) + 0.001/2;
            plot(tt,yy, 'Color', col(ii,:), 'LineWidth', jj);
        end
    end
    ylabel('Firing rate (Hz)')
    xlabel('Time (s)')
    
    if k==1
        title('Inactive session')
    else
        title('Active session')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUA PSTH POOLED

ds = {find(l.coefVar>1.2), find(l.coefVar<0.6)};

for kk=1:2    
    muas = [];
    muas_smooth = [];
    
    for k = 1:length(ds{kk})
        X = l.datasets{ds{kk}(k)}(:,:,:,l.time>-0.05 & l.time<0.15,:) / 0.01;
        X = mean(X,1);
        X = nanmean(X,5);
        mua = squeeze(X);
        muas = cat(4, muas, mua);
        
        for ii = 1:12
            for jj = 1:3
                ntrial = find(~isnan(l.datasets{ds{kk}(k)}(1,ii,jj,1,:)), 1, 'last');
                yy = l.spikes{ds{kk}(k)}(:,ii,jj,:,1:ntrial);
                yy = yy(:);
                yy = yy(~isnan(yy));
                yy = histcounts(yy, [-0.3:0.001:0.3]);
                yy = yy / size(l.spikes{ds{kk}(k)},1) / ntrial;
                gaussKernel = 1/sqrt(2*pi)/5 * exp(-[-30:30].^2/5^2/2);
                yy = conv(yy, gaussKernel, 'same') * 1000;
                yy = yy(251:450);
                
                muas_smooth(ii,jj,:,k) = yy;
            end
        end
    end
    
    mua_smooth = mean(muas_smooth,4);
    
    mua = squeeze(mean(muas,4));
    baseline = mua(:,:,1:5);
    baseline = mean(baseline(:));
    
    subplot(2,4,(kk-1)*4+2)
    hold on
    axis([-0.05 0.15 0 5])
    plot(xlim, baseline*[1 1], 'Color', [.5 .5 .5])
    plot([0 0], ylim, 'Color', [.5 .5 .5])
    plot([0.02 0.02], ylim, '--', 'Color', [.5 .5 .5])
    plot([0.04 0.04], ylim, '--', 'Color', [.5 .5 .5])
    t = l.time(l.time>-0.05 & l.time<0.15);
    
    col = parula12();
    for ii=1:12
        for jj=1:3
%            y = squeeze(mua(ii,jj,:));
            y = squeeze(mua_smooth(ii,jj,:));
            tt = (-0.05:0.001:0.15-.001) + 0.001/2;
            plot(tt, y, 'Color', col(ii,:), 'LineWidth', jj);
        end
    end
    ylabel('Firing rate (Hz)')
    xlabel('Time (s)')
    if kk==1
        title('Inactive sessions')
    else
        title('Active sessions')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % tuning curves

    muaOnset = mean(mua(:,:,8:9),3);
    mua = mean(mua(:,:,6:end),3);
    
    subplot(2,4,(kk-1)*4+3)
    hold on
    
    col = parula12();
    xlim([0 42])
    plot(xlim, baseline*[1 1], 'Color', [.5 .5 .5])
    for j=1:3
        x = (1:12) + (j-1)*14;
        for i=1:12
            plot(x(i), mua(i,j), 'o', 'MarkerSize', 5+j*2, 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'w')
            plot(x(i), muaOnset(i,j), '.', 'Color', col(i,:)*.6+[1 1 1]*.4, 'MarkerSize', 8)
        end
    end
    set(gca,'XTick',[])
    ylim([0 5])
    ylabel('Firing rate (Hz)')
    xlabel('Stimulus type')
    title('0--150ms')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACROSS SESSIONS

for k=1:length(l.datasets)
    %X = l.datasets{k}(:,:,:,l.time>0.02 & l.time<0.04,:) / 0.01;
    X = l.datasets{k}(:,:,:,l.time>0 & l.time<0.15,:) / 0.01;
    X = mean(X,1);
    X = mean(X,4);
    X = nanmean(X,5);
    mua = squeeze(X);
    tuningStrength(k) = mua(1,3)/mua(12,1);
    tuningStrength(k) = min(tuningStrength(k), 5);
end

subplot(2,4,[4 8])
axis([0 2.2 0 5])
axis square
hold on
myscatter(l, tuningStrength)
xlabel('CV')
ylabel('Contra loud / ipsi silent')
[r,p] = corr(l.coefVar(:), tuningStrength(:));
text(0.2, 4.5, ['r=' num2str(r,2) ', p=' num2str(p,3)])

subplot(241)
text(-0.2, 1.15, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(242)
text(-0.2, 1.15, 'B', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(243)
text(-0.2, 1.15, 'C', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(245)
text(-0.2, 1.15, 'D', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(246)
text(-0.2, 1.15, 'E', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(247)
text(-0.2, 1.15, 'F', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
subplot(2,4,[4 8])
text(-0.2, 1.15, 'G', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)


%%%%%%%%%%%%%%%
% export

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureMUA.pdf','-dpdf','-r0')

end
