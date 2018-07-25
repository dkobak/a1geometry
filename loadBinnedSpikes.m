function [Xtrial, Xspikes, shanks, inactLevel, fano, cv] = loadBinnedSpikes(fileName, binBorders, intervalsToUse)

% Loads a file with fileName, bins spikes in between the binBorders
% Returns Xtrial(N, ild, abl, timebin, trial) array

if nargin<3
    intervalsToUse = [-inf inf];
end

datadir = 'data/';
mydata = load([datadir fileName]);

stimInterval = [1 length(mydata.stimOnsetTime)];
% the initial stimuli in the 1st file are some other stuff
if length(mydata.stimOnsetTime) == 3648
    stimInterval(1) = 49;
end

[t,f,fano,cv] = getActivationCurve(mydata);
ibeg = find(t>mydata.stimOnsetTime(stimInterval(1)),1);
iend = find(t>mydata.stimOffsetTime(stimInterval(2)),1)-1;
if isempty(iend)
    iend = length(t);
end
for i = ibeg : iend
    if ~isempty(find(sum(intervalsToUse-t(i)/60 > 0, 2) == 1, 1))
        tin(i) = 1;
    end
end
inactLevel = mean(f(tin==1));
fano = mean(fano(tin==1));
cv = mean(cv(tin==1));

shanks = mydata.shankIndex_fNeuron;
        
[~,~,ablG] = unique(mydata.ABL);
[~,~,ildG] = unique(mydata.ILD);

XevokedT = nan(length(mydata.spikes), max(ildG), max(ablG), length(binBorders)-1, 200);
nn = zeros(max(ildG), max(ablG));
Xspikes = nan(length(mydata.spikes), max(ildG), max(ablG), 100, 200);
maxnspikes = 0;

for i = stimInterval(1):stimInterval(2)
    if ~isempty(find(sum(intervalsToUse-mydata.stimOnsetTime(i)/60 > 0, 2) == 1, 1))
        nn(ildG(i), ablG(i)) = nn(ildG(i), ablG(i)) + 1;
        for n = 1:length(mydata.spikes)
            h = histc(mydata.spikes{n}, mydata.stimOnsetTime(i) + binBorders);
            XevokedT(n, ildG(i), ablG(i), :, nn(ildG(i), ablG(i))) = h(1:end-1);
            
            hh = mydata.spikes{n}-mydata.stimOnsetTime(i);
            hh = hh((hh>binBorders(1)-0.15) & hh<binBorders(end)+0.15);
            Xspikes(n, ildG(i), ablG(i),  1:length(hh), nn(ildG(i), ablG(i))) = hh;
            if length(hh) > maxnspikes
                maxnspikes = length(hh);
            end
        end
    end
end

Xtrial = XevokedT(:,:,:,:,1:max(nn(:)));
Xspikes = Xspikes(:,:,:,1:maxnspikes,1:max(nn(:)));

