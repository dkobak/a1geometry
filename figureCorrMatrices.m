function figureCorrMatrices(l, covorcorr)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

% cov or corr
if nargin==2 && strcmp(covorcorr, 'cov')
    covorcorr = 'cov';
    colorrange = {[-2 2], [-1 1]};
    ylims = [.3 .8 15
             .3 .8 15];    
    ylimspec = {[0 .6], [0 .25]};
else
    covorcorr = 'corr';
    colorrange = {[-1 1], [-.25 .25]};
    ylims = [.3 .5 10
             .3 .5 40];
    ylimspec = {[0 .3], [0 .1]};
end

fig = figure('Position', [100 100 1500 500]);
colormap(parula)

numRep = 100; % number of shuffles

%%%%%%%%%%% TWO EXAMPLE DATASETS

examples = [1 7];

linecolors = {[208,28,139]/256, [77,172,38]/256};

for ex = 1:2
    datasetNum = examples(ex);
    
    xx = l.datasets{datasetNum}(:,:,:,l.time>0 & l.time<0.15,:);
    xx = squeeze(sum(xx,4)) / 0.15;
    Xsignal = nanmean(xx,4);
    Xnoise = bsxfun(@minus, xx, Xsignal);

    % subtracting running mean
    window = 11;
    for a=1:size(xx,2)
        for b=1:size(xx,3)
            n = sum(~isnan(xx(1,a,b,:)));
            for w=1:n
                win = [max(1, w-(window-1)/2) min(n, w+(window-1)/2)];
                Xnoise(:,a,b,w) = bsxfun(@minus, xx(:,a,b,w), mean(xx(:,a,b,win),4));
            end
        end
    end    
    
    Xsignal = Xsignal(:,:)';
    Xnoise = Xnoise(:,:)';
    Xnoise = Xnoise(~isnan(Xnoise(:,1)), :);
    Xall = xx(:,:)';
    Xall = Xall(~isnan(Xall(:,1)),:);

    v = var(Xnoise) + var(Xsignal);
    Xsignal = Xsignal(:, v>0);
    Xnoise  = Xnoise(:, v>0);
    %Xall = Xall(:, v>0);
    
    Xsignal = bsxfun(@minus, Xsignal, mean(Xsignal));
    Xnoise =  bsxfun(@minus, Xnoise,  mean(Xnoise));
    if strcmp(covorcorr, 'corr')
        Xsignal = bsxfun(@times, Xsignal, 1./(std(Xsignal)));
        Xnoise =  bsxfun(@times, Xnoise,  1./(std(Xnoise) ));
    end
    
    X = {Xsignal, Xnoise};
    Csignal = cov(Xsignal);
    Cnoise  = cov(Xnoise);
    C = {Csignal, Cnoise};
        
    for mode = 1:2 % signal and noise
        [v, e] = eig(C{mode});
        spectrum{ex,mode} = sort(abs(diag(e)), 'descend') / sum(abs(diag(e)));
        
        % using the same sorting for Csignal and Cnoise (based on beta_ild)
        if mode == 1 
            beta = getBetas(l, datasetNum);
            [~,neuronOrder{ex}] = sort(beta(:,1), 'ascend');
        end
        
        subplot(2,5,(ex-1)*1 + (mode-1)*5 + 1)
        imagesc(C{mode}(neuronOrder{ex},neuronOrder{ex}), colorrange{mode})
        axis square
        set(gca, 'YTick', [], 'XTick', [])        
        
        esh{ex,mode} = zeros(numRep, size(X{mode},2));
        for rep = 1:numRep
            xsh = X{mode};
            for i=1:size(xsh,2)
                xsh(:,i) = xsh(randperm(size(xsh,1)),i);
            end
            esh{ex,mode}(rep,:) = sort(abs(eig(cov(xsh))), 'descend');
            esh{ex,mode}(rep,:) = esh{ex,mode}(rep,:)/sum(esh{ex,mode}(rep,:));
        end
        
        dimens(ex, mode) = 0;
    end
end

subplot(252)
title('Active session')
subplot(251)
title('Inactive session')
ylabel('Signal', 'fontweight', 'bold', 'fontsize', 12)
subplot(257)
title('Active session')
subplot(256)
title('Inactive session')
ylabel('Noise', 'fontweight', 'bold', 'fontsize', 12)


%%%%%%%%% ACROSS SESSIONS

for datasetNum = 1:length(l.datasets)
    fprintf('.')
   
    xx = l.datasets{datasetNum}(:,:,:,l.time>0 & l.time<0.15,:);
    xx = squeeze(sum(xx,4)) / 0.15;
    Xsignal = nanmean(xx,4);
    Xnoise = bsxfun(@minus, xx, Xsignal);
    
    % subtracting running mean
    window = 11;
    for a=1:size(xx,2)
        for b=1:size(xx,3)
            n = sum(~isnan(xx(1,a,b,:)));
            for w=1:n
                win = [max(1, w-(window-1)/2) min(n, w+(window-1)/2)];
                Xnoise(:,a,b,w) = bsxfun(@minus, xx(:,a,b,w), mean(xx(:,a,b,win),4));
            end
        end
    end   
    
    Xsignal = Xsignal(:,:)';
    Xnoise = Xnoise(:,:)';
    Xnoise = Xnoise(~isnan(Xnoise(:,1)), :);
    
    v = var(Xnoise) + var(Xsignal);
    Xsignal = Xsignal(:, v>0);
    Xnoise  = Xnoise(:, v>0);
    
    Xsignal = bsxfun(@minus, Xsignal, mean(Xsignal));
    Xnoise =  bsxfun(@minus, Xnoise,  mean(Xnoise));
    if strcmp(covorcorr, 'corr')
        Xsignal = bsxfun(@times, Xsignal, 1./(std(Xsignal)));
        Xnoise =  bsxfun(@times, Xnoise,  1./(std(Xnoise)));
    end
    
    X = {Xsignal, Xnoise};
    
    Csignal = cov(Xsignal);
    Cnoise  = cov(Xnoise);
    C = {Csignal, Cnoise};
    
    for mode = 1:2 % signal and noise
        spectrum = sort(abs(eig(C{mode})), 'descend');
        toppc(datasetNum, mode) = spectrum(1)/sum(spectrum);
        meancorr(datasetNum, mode) = mean(C{mode}(eye(size(C{mode},1))==0)) / ...
                                     mean(C{mode}(eye(size(C{mode},1))==1));
        
        esh = zeros(numRep, size(X{mode},2));
        for rep = 1:numRep
            xsh = X{mode};
            for i=1:size(xsh,2)
                xsh(:,i) = xsh(randperm(size(xsh,1)),i);
            end
            esh(rep,:) = sort(abs(eig(cov(xsh))), 'descend');
        end

        numdim(datasetNum, mode) = find(spectrum' < prctile(esh,99), 1) - 1;
    end
end
fprintf('\n')

titles = {'Average correlation', 'Fraction of var in PC1', 'Dimensionality'};
data = {meancorr, toppc, numdim};
textverposition = [.8 .1 .8; .8 .8 .1];
for mode = 1:2
    for stat = 1:3
        subplot(2,5,(mode-1)*5 + 2 + stat)
        
        y = data{stat}(:,mode);
        myscatter(l, y);

        ylim([0 ylims(mode,stat)])
        ylabel(titles{stat})
        xlabel('CV')
        
        [r,p] = corr(l.coefVar(:), y(:));
        if (p>0.04 && p<0.06)
            precision = 2;
        else
            if (p<0.0001)
                precision = '%.5f';
            else
                precision = 1;
            end
        end
        text(.2, ylims(mode,stat)*textverposition(mode,stat), ...
            ['$$r=' num2str(r,2) ',\, p=' num2str(p,precision) '$$'], 'Interpreter', 'latex')
    end
end

letters = 'ABCDEFGHIJ';
xd = [-.15 -.15 -.3 -.3 -.3 -.15 -.15 -.3 -.3 -.3];
for i = 1:10
    subplot(2,5,i)
    text(xd(i), 1.12, letters(i), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

axes('Position', [0.24 0.65 0.01 0.2])
set(gca, 'visible', 'off')
caxis(colorrange{1})
colorbar
axes('Position', [0.24 0.18 0.01 0.2])
set(gca, 'visible', 'off')
caxis(colorrange{2})
colorbar

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'figures/figureCorrMatrices.pdf','-dpdf','-r0')

end

% This is needed to order the neurons when displaying corr matrices
function beta = getBetas(l,d)
    X = squeeze(sum(l.datasets{d}(:,:,:,l.time>0 & l.time<0.15,:),4) / 0.15);
    ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
    abl = [20 40 60];
    abl = abl-40;

    beta = [];
    for n = 1:size(X,1)
        fr = X(n,:,:,:);
        ind = ~isnan(fr(:));
        ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
        ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
        b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)]);
        beta = [beta; b'];
    end
end

