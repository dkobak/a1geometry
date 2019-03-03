function figureConditioning(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

figure('Position', [100 100 1500 750])

sessions = find(l.coefVar>1.2);

% conditioning on UP or on DOWN
for mode = 1:2       
    muas = [];
    muas_smooth = zeros(size(l.datasets{1},2), size(l.datasets{1},3), 200, length(sessions));
    betas = [];
    
    % go over all inactive datasets
    for k = 1:length(sessions)
        Xbefore = l.datasets{sessions(k)}(:,:,:,l.time>-0.05 & l.time<0,:);
        Xbefore = sum(Xbefore,4);
        Xbefore = sum(Xbefore,1);
        
        X = l.datasets{sessions(k)}(:,:,:,l.time>-0.05 & l.time<0.15,:) / 0.01;
        
        for i=1:size(X,2)
            for j=1:size(X,3)
                if mode == 1
                    X(:,i,j,:, Xbefore(1,i,j,1,:) == 0) = nan;
                else
                    X(:,i,j,:, Xbefore(1,i,j,1,:) > 0) = nan;
                end
                
                % smoothing
                ind = find(~isnan(X(1,i,j,1,:)));
                yy = l.spikes{sessions(k)}(:,i,j,:,ind);
                yy = yy(:);
                yy = yy(~isnan(yy));
                yy = histcounts(yy, [-0.3:0.001:0.3]);
                yy = yy / size(l.spikes{sessions(k)},1) / length(ind);
                gaussKernel = 1/sqrt(2*pi)/5 * exp(-[-30:30].^2/5^2/2);
                yy = conv(yy, gaussKernel, 'same') * 1000;
                yy = yy(251:450);
                muas_smooth(i,j,:,k) = yy;
            end
        end
        
        % getting betas for each neuron
        for n = 1:size(X,1)
            fr = squeeze(sum(X(n,:,:,6:end,:), 4));
            fr = nanmean(fr, 3); % regression on the tuning curves
            ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
            abl = [20 40 60];
            abl = abl-40;
            ildAll = bsxfun(@times, ild', ones([1 3 size(fr,3)]));
            ablAll = bsxfun(@times, abl,  ones([12 1 size(fr,3)]));
            ind = ~isnan(fr(:));
            [b,~,~,~,stats] = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)]);
            if stats(3)<0.01
                betas = [betas; b(1:2)'];
            end
        end
               
        X = mean(X,1);
        X = nanmean(X,5);
        mua = squeeze(X);
        muas = cat(4, muas, mua);        
    end
    
    mua_smooth = mean(muas_smooth, 4);
        
    % plotting mean PSTHs
    subplot(2,3,(mode-1)*3+1)
    hold on
    col = parula12();
    baseline = mua(:,:,1:5);
    baseline = mean(baseline(:));
    axis([-0.05 0.15 0 5])
    plot(xlim, baseline*[1 1], 'Color', [.5 .5 .5])
    plot([0 0], ylim, 'Color', [.5 .5 .5])
    for ii=1:12
        for jj=1:3
            y = squeeze(mua_smooth(ii,jj,:));
            tt = (-0.05:0.001:0.15-.001) + 0.001/2;
            plot(tt, y, 'Color', col(ii,:), 'LineWidth', jj);
        end
    end
    ylabel('Firing rate (Hz)')
    xlabel('Time (s)')
    ylim([0, 5])
    %yticks([0, .01, .02, .03])
    if mode==1
        title('Inactive up-state')
    else
        title('Inactive down-state')
    end
    
    % plotting tuning curves
    subplot(2,3,(mode-1)*3+2)
    hold on

    mua = mean(muas,4);
    baseline = mua(:,:,1:5);
    baseline = mean(baseline(:));
    muaOnset = mean(mua(:,:,8:9),3);
    mua = mean(mua(:,:,6:end),3);

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
    
    % plotting betas
    subplot(2,3,(mode-1)*3+3)
    axis([-1 1 -1 1]*0.1)
    axis square
    hold on
    plot([0 0],ylim, 'k')
    plot(xlim,[0 0], 'k')
    xlabel('$\beta_\mathrm{ILD}$','Interpreter','latex')
    ylabel('$\beta_\mathrm{ABL}$','Interpreter','latex')
    scatter(betas(:,1), betas(:,2), 2, [217,95,2]/256, 'MarkerFaceColor', [217,95,2]/256);
    set(gca, 'XTick', [-.1 -0.05 0 0.05 .1], 'YTick', [-.1 -0.05 0 0.05 .1])
end

letters = 'ABCDEF';
for i = 1:6
    subplot(2,3,i)
    text(-0.25, 1.1, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

%%%%%%%%%%%%%%%
% export

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureConditioning.pdf','-dpdf','-r0')

end
