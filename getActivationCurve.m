function [t, f, fano, cv] = getActivationCurve(mydata)

allspikes = sort(vertcat(mydata.spikes{:}));

firstspike = allspikes(1);
lastspike  = allspikes(end);

bigbin  = 30;
minibin = 0.02;

allspikes = sort(vertcat(mydata.spikes{:}));
sumcount = histc(allspikes, 1:minibin:allspikes(end));

t = 1:minibin:allspikes(end);
tt = 1:bigbin:allspikes(end);
numOfSilentBins = [];
for ttt = 1:length(tt)-1
    x = sumcount(t>=tt(ttt) & t<tt(ttt+1));

    numOfSilentBins(ttt) = sum(x==0);
    fano(ttt) = var(x)/mean(x);
    cv(ttt) = std(x)/mean(x);
end

t = tt(1:end-1);
f = numOfSilentBins/(bigbin/minibin)*100;
