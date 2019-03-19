function preprocessSpikecounts()

time = [-0.025 0.025 0.0775 0.125];

outputFileName = 'evokedResponses.mat';
datadir = 'spikecounts/';

files = dir(datadir);
fileList = {};
for i=3:length(files)    % note that first two are '.' and '..'
    if strcmp(files(i).name(end-3:end), '.csv')
        fileList(length(fileList)+1) = {files(i).name};
    end
end

for filenum = 1:length(fileList)
    display(['Processing file ' num2str(filenum) '...'])
    
    X = csvread([datadir fileList{filenum}]);
    
    for u = unique(X(2:end,1))'
        rows = X(:,1) == u;
        Xpart = X(rows,:);
        
        num = Xpart(1, 1);
        shanks{num} = X(1, 7:end);
        coefVar(num,1) = Xpart(1, 2);
        
        [~,~,ablG] = unique(Xpart(1:4:end, 4));
        [~,~,ildG] = unique(Xpart(1:4:end, 5));
        Xtrial = nan(size(X,2)-6, max(ildG), max(ablG), 4, 200);
        nn = zeros(max(ildG), max(ablG));
        for i = 1:size(Xpart,1)/4
            nn(ildG(i), ablG(i)) = nn(ildG(i), ablG(i)) + 1;
            Xtrial(:, ildG(i), ablG(i), :, nn(ildG(i), ablG(i))) = Xpart((i-1)*4 + 1 : (i-1)*4 + 4, 7:end)';
        end
        Xtrial = Xtrial(:,:,:,:,1:max(nn(:)));
        datasets{num} = Xtrial;
    end

    save(outputFileName, 'datasets', 'time', 'shanks', 'coefVar')
end
