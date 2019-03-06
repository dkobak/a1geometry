function summaryActivations()

datadir = 'data/';

files = dir(datadir);
fileList = {};
for i=3:length(files)    % note that first two are '.' and '..'
    if strcmp(files(i).name(end-3:end), '.mat')
        fileList(length(fileList)+1) = {files(i).name};
    end
end

fig = figure('Position', [100 100 1500 900]);

for fileNum = 1:length(fileList)
    pause(0.01)
    
    mydata = load([datadir fileList{fileNum}]);
    
%     % removing one weird neuron that fires all the time
%     if strcmp(fileList{fileNum}, 'AC02-EXP21_Block-1_dataSet.mat')
%         mydata.spikes = mydata.spikes([1:16 18:end]);
%     end
    
    [t, f, fano, cv] = getActivationCurve(mydata);
    
    subplot(4,5,fileNum)
    u = find(fileList{fileNum}=='_', 1, 'last');
    title(fileList{fileNum}(6:u-1),'interpreter','none')
    set(gca, 'FontName', 'Ubuntu')
    hold on
    plot(t/60, f)
    axis([0 t(end)/60+1 0 90])
    plot(mydata.stimOnsetTime(1) * [1 1]/60, ylim, 'k')
    plot(mydata.stimOnsetTime(end) * [1 1]/60, ylim, 'k')
    text(5,70,['Overall mean: ' num2str(mean(f),3) '%'])
    
    plot(t/60, cv*40, 'k')
end
