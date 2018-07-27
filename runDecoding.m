function runDecoding(l)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

path(path, 'glmnet_matlab/')

for d = 1:length(l.datasets)
    fprintf(['Dataset #' num2str(d) '\n'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DECODING ONLY CONTRA - LINEAR
    fprintf('linear decoding of contra: ')
    for abl = 1:3
        X = squeeze(sum(l.datasets{d}(:, 1:6, abl, l.time>0 & l.time<.15, :), 4));
        y = bsxfun(@times, ones(size(X)), [20 10 6 4.5 3 1.5]);
        [X,y] = cleanup(X,y);
        decoding_abl_contra_linear(d,abl) = nestedcv(X,y,[],'linear');
        fprintf('*')
    end
    fprintf('\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DECODING ONLY IPSI - LINEAR
    fprintf('linear decoding of ipsi: ')
    for abl = 1:3
        X = squeeze(sum(l.datasets{d}(:, 7:12, abl, l.time>0 & l.time<.15, :), 4));
        y = bsxfun(@times, ones(size(X)), [1.5 3 4.5 6 10 20]);
        [X,y] = cleanup(X,y);
        decoding_abl_ipsi_linear(d,abl) = nestedcv(X,y,[],'linear');
        fprintf('*')
    end
    fprintf('\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMON DECODER FOR SIGN(ILD) 
    fprintf('logistic regression for sign(ILD): ')
    X = squeeze(sum(l.datasets{d}(:, :, :, l.time>0 & l.time<.15, :), 4));
    y = bsxfun(@times, ones([1 12 3 size(X,4)]), [-1 -1 -1 -1 -1 -1 1 1 1 1 1 1]);
    group = repmat([(1:12)' (13:24)' (25:36)'], [1 1 size(X,4)]);
    [X,y,group] = cleanup(X,y,group);
    a = nestedcv(X,y,group);
    decoding_psych(d,1:12,1:3) = reshape(a, [12 3]);
    fprintf('\n')

    % Weights

    myoptions = [];
    myoptions.alpha = 0.5;
    cverr = cvglmnet(X, y>0, 'binomial', myoptions);
    decoderWeights{d} = cverr.glmnet_fit.beta(:, cverr.lambda==cverr.lambda_min);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMON DECODER FOR SIGN(ILD) - TIME RESOLVED, ACT ONLY
    if l.coefVar(d) < 0.6
        fprintf('Active session. Performing time-resolved decoding (15 steps): ')
        for t = 1:15
            X = squeeze(sum(l.datasets{d}(:, :, :, l.time>0 & l.time<t/100, :), 4));
            y = bsxfun(@times, ones([1 12 3 size(X,4)]), [-1 -1 -1 -1 -1 -1 1 1 1 1 1 1]);
            group = repmat([(1:12)' (13:24)' (25:36)'], [1 1 size(X,4)]);
            [X,y,group] = cleanup(X,y,group);
            a = nestedcv(X,y,group);
            decoding_psych_intime(d, 1:12, 1:3, t) = reshape(a, [12 3]);
            fprintf('*')
        end
        fprintf('\n')
    else
        decoding_psych_intime(d, 1:12, 1:3, 1:15) = nan;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ABL-dependent DECODER FOR SIGN(ILD) 
    fprintf('logistic regression for sign(ILD), ABL dependent: ')
    for abl = 1:3
        X = squeeze(sum(l.datasets{d}(:, :, abl, l.time>0 & l.time<.15, :), 4));
        y = bsxfun(@times, ones([1 12 size(X,3)]), [-1 -1 -1 -1 -1 -1 1 1 1 1 1 1]);
        group = repmat((1:12)', [1 size(X,3)]);
        [X,y,group] = cleanup(X,y,group);
        a = nestedcv(X,y,group);
        decoding_psych_perabl(d,1:12,abl) = a;
        fprintf('*')
    end
    fprintf('\n')
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMON DECODER FOR SIGN(ILD) WITH SHUFFLE
    fprintf('logistic regression for sign(ILD), shuffled (10 reps): ')
    for rep = 1:10
        datasetSh = shuffle(l.datasets{d}, rep);
        X = squeeze(sum(datasetSh(:, :, :, l.time>0 & l.time<.15, :), 4));
        y = bsxfun(@times, ones([1 12 3 size(X,4)]), [-1 -1 -1 -1 -1 -1 1 1 1 1 1 1]);
        group = repmat([(1:12)' (13:24)' (25:36)'], [1 1 size(X,4)]);
        [X,y,group] = cleanup(X,y,group);
        a = nestedcv(X,y,group);
        decoding_psych_shuffle(d,1:12,1:3,rep) = reshape(a, [12 3]);
        fprintf('*')
    end
    fprintf('\n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMON DECODER FOR SIGN(ILD) WITH SHUFFLE ONLY TRAINING
    fprintf('logistic regression for sign(ILD), shuffled training (10 reps): ')
    for rep = 1:10
        X = squeeze(sum(l.datasets{d}(:, :, :, l.time>0 & l.time<.15, :), 4));
        y = bsxfun(@times, ones([1 12 3 size(X,4)]), [-1 -1 -1 -1 -1 -1 1 1 1 1 1 1]);
        group = repmat([(1:12)' (13:24)' (25:36)'], [1 1 size(X,4)]);
        [X,y,group] = cleanup(X,y,group);
        a = nestedcv(X,y,group,[],'shuffle_training',rep);
        decoding_psych_shuffleTraining(d,1:12,1:3,rep) = reshape(a, [12 3]);
        fprintf('*')
    end
    fprintf('\n')
        
    save('decoding_results.mat', 'decoding_abl_contra_linear', 'decoding_abl_ipsi_linear', ...
         'decoding_psych', 'decoding_psych_perabl', 'decoding_psych_shuffle', ...
         'decoding_psych_shuffleTraining', 'decoderWeights', 'decoding_psych_intime')
end

end

function [XX, yy, gg] = cleanup(X,y,group)
    XX = X(:,:)';
    yy = y(:);
    indx = ~isnan(XX(:,1));
    XX = XX(indx,:);
    yy = yy(indx);
    XX = bsxfun(@minus, XX, mean(XX));
    
    if nargin==3
        gg = group(:);
        gg = gg(indx);
    end
end

function dataSh = shuffle(data, rep)
    rng(rep)
    
    dataSh = nan(size(data));
    for j=1:size(data,2)
        for k=1:size(data,3)
            numTr = find(~isnan(squeeze(data(1,j,k,1,:))), 1, 'last');
            for i=1:size(data,1)
                permInd = randperm(numTr);
                dataSh(i,j,k,:,1:numTr) = data(i,j,k,:,permInd);
            end
        end
    end
end

function accuracy = nestedcv(X, y, group, iflinear, ifshuffle, shuffleseed)
    % for glmnet
    myoptions = [];
    myoptions.alpha = 0.5;
    
    % `group` argument allows to measure test performance on several groups
    if nargin < 3 || isempty(group)
        group = ones(size(y));
    end
    
    if nargin < 4 || isempty(iflinear)
        iflinear = 'logistic';
    end
    
    if nargin < 5 || isempty(ifshuffle)
        ifshuffle = 'none';
    end
    
    rng(42)
    ind = randperm(size(X,1));
    X = X(ind,:);
    y = y(ind);
    group = group(ind);
    foldn = floor(size(X,1)/10);
    
    yhatall = y*0;
    
    for fold = 1:10
        fprintf('.')
        
        if fold==10
            endf = size(X,1);
        else
            endf = fold*foldn;
        end
        indtest = (fold-1)*foldn + 1 : endf;
        indtrain = setdiff(1:size(X,1), indtest);
        
        Xtrain = X(indtrain,:);
        ytrain = y(indtrain);
        Xtest  = X(indtest,:);
        ytest  = y(indtest);
        grouptest = group(indtest);
        
        if strcmp(ifshuffle, 'shuffle_training')
            rng(shuffleseed)
            for gg = 1:max(group(indtrain))
                ii = find(group(indtrain) == gg);
                for jj = 1:size(Xtrain,2)
                    Xtrain(ii,jj) = Xtrain(ii(randperm(length(ii))), jj);
                end
            end
        end
        
        if strcmp(iflinear, 'logistic')
            cverr = cvglmnet(Xtrain, ytrain>0, 'binomial', myoptions);
            yhat = cvglmnetPredict(cverr, Xtest, 'lambda_min');
            yhatall(indtest) = yhat;
        else
            cverr = cvglmnet(Xtrain, ytrain, [], myoptions);
            yhat = cvglmnetPredict(cverr, Xtest, 'lambda_min');
            yhatall(indtest) = yhat;
        end
    end
    
    if strcmp(iflinear, 'logistic')
        for g = 1:max(group)
            accuracy(g) = mean(sign(yhatall(group==g))==sign(y(group==g)));
        end
    else
        for g = 1:max(group)
            accuracy(g) = 1 - sum((y(group==g) - yhatall(group==g)).^2) / ...
                              sum((y(group==g) - mean(y(group==g))).^2);
        end
    end
end

% function accuracy = simplecv(X,y)
%     % for glmnet
%     myoptions = [];
%     myoptions.alpha = 0.5;
%     
%     rng(42)
%     cverr = cvglmnet(X, y>0, 'binomial', myoptions);
%     yhat = cvglmnetPredict(cverr, X, 'lambda_min');
%     accuracy = sum(sign(yhat)==sign(y)) / length(y);
% end
