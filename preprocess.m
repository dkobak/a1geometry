function preprocess()

% sessions that need to be split into two

filesToSplit{1} = 'AC02-EXP08_Block-1_dataSet.mat';
filesToSplit{2} = 'AC02-EXP09_Block-1b_dataSet.mat';
filesToSplit{3} = 'AC02-EXP23_Block-2_dataSet.mat';

actPeriods{1} = [76 100];
actPeriods{2} = [0 10; 12 15; 46 58];
actPeriods{3} = [0 26];

inactPeriods{1} = [0 76];
inactPeriods{2} = [10 12; 15 46; 58 70];
inactPeriods{3} = [26 70];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datasets = {};
inactLevel = [];
fanoFactor = [];
coefVar = [];
ifSplit = [];
fileNum = [];
fileName = {};
expNum = [];
blockNum = [];

timeEdge = -0.15:0.01:0.15;
time = timeEdge(1:end-1) + mean(diff(timeEdge))/2;

outputFileName = 'evokedResponses_150to150.mat';

datadir = 'data/';

files = dir(datadir);
fileList = {};
for i=3:length(files)    % note that first two are '.' and '..'
    if strcmp(files(i).name(end-3:end), '.mat')
        fileList(length(fileList)+1) = {files(i).name};
    end
end

tic
for filenum = 1:length(fileList)
    display(['Processing file ' num2str(filenum) '...'])
    
    i = find(strcmp(filesToSplit, fileList{filenum}));
    if isempty(i)
        [Xtrial, Xspikes, sh, il, fano, cv, Xspikecount] = loadBinnedSpikes(fileList{filenum}, timeEdge);
        fileNum = [fileNum; filenum];
        fileName{length(fileNum)} = fileList{filenum};
        ifSplit = [ifSplit; 0];
        expNum = [expNum; str2num(fileList{filenum}(9:10))];
        blockNum = [blockNum; str2num(fileList{filenum}(18))];
        inactLevel = [inactLevel; il];
        fanoFactor = [fanoFactor; fano];
        coefVar = [coefVar; cv];
        datasets{length(fileNum)} = Xtrial;
        spikes{length(fileNum)} = Xspikes;
        shanks{length(fileNum)} = sh;
        
        Xspikecount(2:end,1) = length(fileNum);
        Xspikecount(1, 1:6) = 0;
        Xspikecount(1, 7:end) = sh;
    else
        display('Splitting in two...')
        
        [Xtrial, Xspikes, sh, il, fano, cv, Xspikecount1] = loadBinnedSpikes(fileList{filenum}, timeEdge, actPeriods{i});
        fileNum = [fileNum; filenum];
        fileName{length(fileNum)} = fileList{filenum};
        ifSplit = [ifSplit; 1];
        expNum = [expNum; str2num(fileList{filenum}(9:10))];
        blockNum = [blockNum; str2num(fileList{filenum}(18))];
        inactLevel = [inactLevel; il];
        fanoFactor = [fanoFactor; fano];
        coefVar = [coefVar; cv];
        datasets{length(fileNum)} = Xtrial;
        spikes{length(fileNum)} = Xspikes;
        shanks{length(fileNum)} = sh;
        
        [Xtrial, Xspikes, sh, il, fano, cv, Xspikecount2] = loadBinnedSpikes(fileList{filenum}, timeEdge, inactPeriods{i});
        fileNum = [fileNum; filenum];
        fileName{length(fileNum)} = fileList{filenum};
        ifSplit = [ifSplit; 1];
        expNum = [expNum; str2num(fileList{filenum}(9:10))];
        blockNum = [blockNum; str2num(fileList{filenum}(18))];
        inactLevel = [inactLevel; il];
        fanoFactor = [fanoFactor; fano];
        coefVar = [coefVar; cv];
        datasets{length(fileNum)} = Xtrial;  
        spikes{length(fileNum)} = Xspikes;
        shanks{length(fileNum)} = sh;
        
        Xspikecount1(2:end,1) = Xspikecount1(2:end,2)*0 + length(fileNum) - 1;
        Xspikecount2(2:end,1) = Xspikecount2(2:end,2)*0 + length(fileNum);
        Xspikecount = nansum(cat(3,Xspikecount1, Xspikecount2),3);
        Xspikecount(1, 1:6) = 0;
        Xspikecount(1, 7:end) = sh;
    end
    
    csvwrite(['spikecounts/' fileList{filenum}(1:end-3) 'csv'], Xspikecount)
    
    save(outputFileName, 'datasets', 'time', ...
        'inactLevel', 'ifSplit', 'fileNum', 'fileName', 'expNum', 'blockNum', ...
        'shanks', 'fanoFactor', 'coefVar', 'spikes', '-v7.3')
end
toc
