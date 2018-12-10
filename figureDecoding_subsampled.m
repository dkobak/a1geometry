function figureDecoding_subsampled(l)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

load decoding_results.mat

colAct = [27,158,119]/256;
colInact = [217,95,2]/256;

ilds = 1:12;

figure('Position', [100 100 1000 600])
letters = 'ABC'; 
for abl=1:3
    subplot(2,3,abl)
    hold on
    avAcc = squeeze(mean(decoding_psych_intime(:,ilds,abl,:),2));
    p1 = plot((1:15)*10, avAcc(l.coefVar<0.6,:), 'Color', colAct);
    p2 = plot((1:15)*10, avAcc(l.coefVar>1.2,:), 'Color', colInact);
    p3 = plot((1:15)*10, avAcc(l.coefVar>0.6 & l.coefVar<1.2,:), 'Color', [.6, .6, .6]);
    p4 = plot((1:15)*10, mean(avAcc(l.coefVar<0.6,:)), 'Color', colAct, 'LineWidth', 3);
    p5 = plot((1:15)*10, mean(avAcc(l.coefVar>1.2,:)), 'Color', colInact, 'LineWidth', 3);
    
    xlabel('Window end (ms)')
    ylabel('Average accuracy')
    xlim([0 160])
    ylim([.45 1])
    title(['ABL = ' num2str(abl*20)])
    text(-0.2, 1.1, letters(abl), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
   
    if abl==1
        legend([p1(1),p2(1),p3(1),p4,p5], {'Active sessions', 'Inactive sessions', ...
            'Interm. sessions', 'Average active', 'Average inactive'}, ...
            'Position', [.14 .85 .14 .1]);
        legend boxoff
        title('')
    end
end

letters = 'DEF'; 
for abl=1:3
    subplot(2,3,abl+3)
    hold on
    avAcc = squeeze(mean(nanmean(decoding_psych_subsampled(:,ilds,abl,:,:),5),2));
    for session = 1:numel(l.datasets)
        n = size(l.datasets{session}, 1);
        ns = 10:10:ceil(n/10)*10;
        ns(end) = n;
        
        if l.coefVar(session)<0.6
            plot(ns, avAcc(session,1:numel(ns)), 'Color', colAct)    
        elseif l.coefVar(session)>1.2
            plot(ns, avAcc(session,1:numel(ns)), 'Color', colInact)    
        else
            plot(ns, avAcc(session,1:numel(ns)), 'Color', [.6, .6, .6])    
        end
    end    
    
    plot((1:10)*10, nanmean(avAcc(l.coefVar<0.6,1:10)), 'Color', colAct, 'LineWidth', 3)
    plot((1:10)*10, nanmean(avAcc(l.coefVar>1.2,1:10)), 'Color', colInact, 'LineWidth', 3)
    
    xlabel('Number of neurons')
    ylabel('Average accuracy')
    xlim([0 150])
    ylim([.45 1])
    text(-0.2, 1.1, letters(abl), 'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 17)
end

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*1.1])
print(h,'figures/figureDecodingSubsampled.pdf','-dpdf','-r0')


