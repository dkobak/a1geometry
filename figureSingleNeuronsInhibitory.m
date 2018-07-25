function figureSingleNeuronsInhibitory(l)

% THIS IS NOT A FIGURE :-)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

datasetInd = find(l.coefVar < .6);
[beta, dataset, neuron] = getBetas(datasetInd, l.time>0 & l.time<0.15);
[~, indcontra] = sort(beta(:,1));
indcontra = indcontra(1:25);
[~, indipsi] = sort(beta(:,1), 'descend');
indipsi = indipsi(1:25);

count1 = 0;
for i = 1:length(indcontra)
    d = dataset(indcontra(i));
    n = neuron(indcontra(i));
    baseline = squeeze(sum(l.datasets{d}(n,:,:,l.time>-.05 & l.time<0,:),4) / 0.05);
    baseline = nanmean(baseline(:));
    response = squeeze(sum(l.datasets{d}(n,1:6,:,l.time>0 & l.time<0.05,:),4) / 0.05);
    response = nanmean(response,3);
    response = mean(response(:));
    if response < baseline
        count1 = count1+1;
    end    
end

count2 = 0;
for i = 1:length(indipsi)
    d = dataset(indipsi(i));
    n = neuron(indipsi(i));
    baseline = squeeze(sum(l.datasets{d}(n,:,:,l.time>-.05 & l.time<0,:),4) / 0.05);
    baseline = nanmean(baseline(:));
    response = squeeze(sum(l.datasets{d}(n,7:12,:,l.time>0 & l.time<0.05,:),4) / 0.05);
    response = nanmean(response,3);
    response = mean(response(:));
    if response < baseline
        count2 = count2+1;
    end    
end

[count1 count2; length(indcontra) length(indipsi)]
[~,p] = fishertest([count1 count2; length(indcontra) length(indipsi)])

    function [beta, dataset, neuron] = getBetas(datasetInd, timeinterval)
        beta = [];
        dataset = [];
        neuron = [];
        
        ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20];
        abl = [20 40 60] - 40;

        for d = 1:length(datasetInd)
            X  = squeeze(sum(l.datasets{datasetInd(d)}(:,:,:,timeinterval,:),4) / 0.15);
            for n = 1:size(X,1)
                fr = X(n,:,:,:); 
                ind = ~isnan(fr(:));
                ildAll = bsxfun(@times, ild', ones([1 3 size(X,4)]));
                ablAll = bsxfun(@times, abl,  ones([12 1 size(X,4)]));
                b = regress(zscore(fr(ind)), [ildAll(ind) ablAll(ind) ones(sum(ind), 1)], 0.01);
                beta = [beta; b(1:2)'];   
                dataset = [dataset; datasetInd(d)];
                neuron = [neuron; n];
            end
        end
    end

    function scatterBetas(beta, col)
        axis([-1 1 -1 1]*0.065)
        axis square
        hold on
        plot([0 0],ylim, 'k')
        plot(xlim,[0 0], 'k')
        xlabel('ILD')
        ylabel('ABL')
        if nargin==1
            col = [0 0 0];
        end
        scatter(beta(:,1), beta(:,2), 10, col, 'MarkerFaceColor', col);
        set(gca, 'XTick', [-0.05 0 0.05], 'YTick', [-0.05 0 0.05])
    end
end