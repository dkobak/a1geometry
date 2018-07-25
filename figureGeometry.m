function figureGeometry(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

addpath('../dpca/matlab')

begtime = 0.00;
endtime = 0.15;

for d = 1:length(l.datasets)
    X = squeeze(sum(l.datasets{d}(:,:,:,l.time>begtime & l.time<endtime ,:),4) / (endtime-begtime));
    Xpsth = nanmean(X,4);
    Xnoise = bsxfun(@minus, X, Xpsth);
    
    XXpsth = Xpsth(:,:)';
    meanPSTH = mean(XXpsth);
    meanPSTH = meanPSTH'/norm(meanPSTH);
    XXnoise = Xnoise(:,:)';
    XXnoise = XXnoise(~isnan(XXnoise(:,1)),:);
    [~,~,Vnoise] = svds(XXnoise, 1);
    
    Xabl = squeeze(nanmean(Xpsth,2))';
    [~,~,Vabl] = svds(Xabl - mean(Xabl), 1);
    Xild = squeeze(nanmean(Xpsth,3))';
    [~,~,Vild] = svds(Xild - mean(Xild), 1);
    
    % dpca start
    trialNum = zeros(size(X,1),12,3);
    for i = 1:12
        for j = 1:3
            trialNum(:,i,j) = find(~isnan(X(1,i,j,:)), 1, 'last');
        end
    end
    Xav = nanmean(X,4);
%     optimalLambda = dpca_optimizeLambda(Xav, X, trialNum, 'simultaneous', true, 'numRep', 10, 'display', 'no');
    optimalLambda = 1e-5; % for speed!
    Cnoise = dpca_getNoiseCovariance(Xav, X, trialNum, 'simultaneous', true);
    [~,V,whichMarg] = dpca(Xav, 5, 'lambda', optimalLambda, 'Cnoise', Cnoise);
    Vild_dpca = V(:, find(whichMarg==1,1));
    Vabl_dpca = V(:, find(whichMarg==2,1));
    % dpca end
        
    VVV = [meanPSTH Vabl Vild Vnoise];
    S(d,1:4,1:4) = abs(VVV'*VVV);
    for k = 1:4
        S(d,k,5) = cos(subspace(VVV(:,k), VVV(:,2:3)));
        S(d,5,k) = S(d,k,5);
        S(d,5,5) = 1;
    end
    
    VVV_dpca = [meanPSTH Vabl_dpca Vild_dpca Vnoise];
    S_dpca(d,1:4,1:4) = abs(VVV_dpca'*VVV_dpca);
    for k = 1:4
        S_dpca(d,k,5) = cos(subspace(VVV_dpca(:,k), VVV_dpca(:,2:3)));
        S_dpca(d,5,k) = S(d,k,5);
        S_dpca(d,5,5) = 1;
    end
end

[r,p] = corr(l.coefVar, acosd(S(:,2,4)));
fprintf(['Angle between ABL and noise: ' num2str(r), ', p = ' num2str(p) '\n']) 
[r,p] = corr(l.coefVar, acosd(S(:,3,4)));
fprintf(['Angle between ILD and noise: ' num2str(r), ', p = ' num2str(p) '\n']) 
[r,p] = corr(l.coefVar, acosd(S(:,2,1)));
fprintf(['Angle between ABL and mean:  ' num2str(r), ', p = ' num2str(p) '\n']) 
[r,p] = corr(l.coefVar, acosd(S(:,3,1)));
fprintf(['Angle between ILD and mean:  ' num2str(r), ', p = ' num2str(p) '\n']) 

figure('Position', [100 100 1600 800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCATTER PLOTS ACROSS SESSIONS

indd = [2 3; 4 5; 1 4; 1 5];
titles = {'ABL / ILD', 'Noise / Signal plane', 'Mean/Noise', 'Mean / Signal plane'};
letters = 'DEFG';
rhos = [];
rhos_dpca = [];
for i=1:4
    % subplot(3,18,18*2+[(i-1)*3+1 : (i-1)*3+3] + 1)
    subplot(3, 12, 24 + [(i-1)*3 + 1 : (i-1)*3 + 3])
    
    axis([0 2.2 0 90])
    hold on
    axis square
    set(gca, 'XTick', [0 1 2], 'YTick', 0:30:90)
    xlabel('CV')
    if i==1
        ylabel('Angle (degrees)')
    end
    y = acosd(S(:, indd(i,1), indd(i,2)));
    myscatter(l, y)
    axis([0 2.2 0 90])
    [r,p] = corr(l.coefVar(:), y(:));
    text(.1, 10, ['r=' num2str(r,2)])
    text(.1, 20, ['p=' num2str(p,3)])
    title(titles{i})
    text(-.15, 1.2, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
    
    rhos(i) = r;
    rhos_dpca(i) = corr(l.coefVar(:), acosd(S_dpca(:, indd(i,1), indd(i,2))));
end

d = rhos_dpca - rhos;
fprintf(['Max PCA/dPCA difference in correlations: ' num2str(max(abs(d))) '\n'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO EXAMPLE PCA PLOTS

example_ind = [1 7]; % Example sessions (one active, one inactive)
col = parula12();
titles = {'Inactive session', 'Active session'};
letters = 'BC';

%subplot(3,3,[1 4])
subplot(3,12, [1 2 3 4 13 14 15 16]);
axis square
text(-.1, 1.12, 'A', 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
axis off

for ii = 1:2
    subplot(3,3,[2+ii-1 5+ii-1])
%    subplot(3, 12, reshape([4; 16] + [(ii-1)*4 + 1 : (ii-1)*4 + 4], 8, []))
    axis square
    hold on
    title(titles{ii})
    xlabel('Signal plane, ILD axis')
    ylabel('Signal plane, dir. orth. to ILD axis')
    if ii==2
        axis([-30 30 -30 30])
    else
        axis([-42 28 -35 35])
    end
    text(-.1, 1.12, letters(ii), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)

    d = example_ind(ii);    
    X = squeeze(sum(l.datasets{d}(:,:,:,l.time>begtime & l.time<endtime ,:),4) / (endtime-begtime));
    
    % mean subtraction
    X = bsxfun(@minus, X, mean(nanmean(X,4),1));
    
    Xpsth = nanmean(X,4);

    XXpsth = Xpsth(:,:)';
    XX      = X(:,:)';
    mu = mean(XXpsth);
    XXpsth = bsxfun(@minus, XXpsth, mu);
    XX     = bsxfun(@minus, XX, mu);
%     [~,~,Vpsth] = svd(XXpsth, 0);
    
    Xabl = squeeze(nanmean(Xpsth,2))';
    [~,~,Vabl] = svd(Xabl - mean(Xabl), 'econ');
    Xild = squeeze(nanmean(Xpsth,3))';
    [~,~,Vild] = svd(Xild - mean(Xild), 'econ');
    V = [Vild(:,1) Vabl(:,1)];
    V(:,2) = V(:,2) - V(:,1)*(V(:,1)'*V(:,2));
    V(:,2) = V(:,2)/norm(V(:,2));
    %V = V(:, [2 1]);
    
    plot(xlim, [0 0], 'k')
    Vablproj = Vabl(:,1)'*V;
    plot(50*[-Vablproj(1) Vablproj(1)], 50*[-Vablproj(2) Vablproj(2)], 'k')
    
    XXpsth = reshape((XXpsth * V)', 2, 12, 3, []);
    XX     = reshape((XX     * V)', 2, 12, 3, []);
    
    % single trials
    for i=1:12
        for j=1:3
            dd = squeeze(XX(:,i,j,:));
            plot(dd(1,:), dd(2,:), '.', 'Color', col(i,:)+([1 1 1]-col(i,:))*0)
        end
    end
    
    % ellipses
    for i=1:12
        for j=1:3
            d = squeeze(XX(:,i,j,:))';
            h = plotEllipse(d);
            set(h, 'LineWidth', 1)
            set(h, 'Color', col(i,:))
        end
    end
    
    % psths
    for i=1:12
        for j=1:3
            plot(XXpsth(1,i,j), XXpsth(2,i,j), 'o', 'MarkerSize', 8+j*3, ...
                'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'w')
        end
    end
end

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureGeometry.pdf','-dpdf','-r0')

end
