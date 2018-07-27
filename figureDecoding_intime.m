function figureDecoding_intime(l)

% Load the data if not provided
if nargin==0 || isempty(l)
    l = load('evokedResponses_150to150.mat');
end

myfun = @(w,x) (w(1)-w(2))./(1+exp(-w(3)*(x-w(4)))) + w(2);
w0 = [1 0 1 0];
ild = [-20 -10 -6 -4.5 -3 -1.5 1.5 3 4.5 6 10 20]';

load decoding_results.mat
decoding_psych_intime(:,1:6,:,:) = 1 - decoding_psych_intime(:,1:6,:,:);
decoding_psych_intime = squeeze(mean(decoding_psych_intime(l.coefVar<0.6, :,:,:)));

warning('error', 'stats:nlinfit:IllConditionedJacobian');
warning('error', 'stats:nlinfit:IterationLimitExceeded');
warning('error', 'stats:nlinfit:ModelConstantWRTParam');
warning('error', 'MATLAB:rankDeficientMatrix');
warnings = zeros(1,15);
for t = 1:15
    for d = 1
        for a = 1:3
            D = decoding_psych_intime(:,a,t);
            try
                w = nlinfit(ild, D(:), myfun, w0);
            catch
                w = [nan nan nan nan];
                warnings(t) = warnings(t) + 1;
            end
            x = -20:.1:20;
            y = myfun(w,x);
            y(x<0) = 1-y(x<0);
            integralPerformances(a,t) = mean(y);
        end
    end
end

figure('Position', [100 100 500 300])
hold on
for a=1:3
    plot((1:15)*10, integralPerformances(a,:), 'k', 'LineWidth', 1+a/2)
end
xlabel('Window end (ms)')
ylabel('Integral accuracy')
%title('Mean over 3 active sessions')
xlim([0 160])
ylim([.5 1])
legend({'ABL = 20 Hz', 'ABL = 40 Hz', 'ABL = 60 Hz'}, 'Location', 'SouthEast')
legend boxoff

h = gcf();
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)*1.1])
print(h,'figures/figureDecodingTimeResolved.pdf','-dpdf','-r0')


