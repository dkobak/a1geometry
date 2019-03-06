function makeMovie(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

begtime = 0.00;
endtime = 0.15;

example_ind = [1 7];
labels = {'Inactive session', 'Active session'};
col = parula12();

figure('Position', [0 0 1000 1000]); 
set(gcf,'Color','k')    
axes('Units', 'normalized', 'Position', [0 0 1 1])
set(gca,'Color','k')    
axis(70*[-1 1 -1 1 -1 1])
set(gcf,'Renderer','zbuffer') 
hold on
axis vis3d
axis off

ann = annotation('textbox',[0.45 0.6 0.1 0.1], 'string', ' ', 'Color', 'w', ...
                 'HorizontalAlignment', 'center');
annotation('ellipse', [.495 .495 .01 .01], 'FaceColor', [.5 .5 .5], 'Color', 'w')
             
vidObj = VideoWriter('animations/geometry.avi');
open(vidObj);

for ii = 1:2
    cla
    
    d = example_ind(ii);    
    X = squeeze(sum(l.datasets{d}(:,:,:,l.time>begtime & l.time<endtime ,:),4) / (endtime-begtime));
    XX = X(:,:)';
    
    Xpsth = nanmean(X,4);
    mu = mean(Xpsth(:,:),2);
    Xabl = squeeze(mean(Xpsth,2))';
    [~,~,Vabl] = svd(Xabl - mean(Xabl), 'econ');
    Xild = squeeze(nanmean(Xpsth,3))';
    [~,~,Vild] = svd(Xild - mean(Xild), 'econ');
    V = gramschmidt([mu Vabl(:,1) -Vild(:,1)]);
    
    XX  = reshape((XX * V)', 3, 12, 3, []);

    % ellipsoids
    for ild=1:12
        for abl=1:3
            data = squeeze(XX(:,ild,abl,:))';
            h = plotEllipsoid(data);
            set(h, 'LineStyle', 'none')
            set(h, 'FaceColor', col(ild,:))
            set(h, 'FaceAlpha', .3)
        end
    end
    
    % psths
    for ild=1:12
        for abl=1:3
            data = squeeze(XX(:,ild,abl,:));
            mu = nanmean(data,2);
            plot3(mu(1), mu(2), mu(3), 'o', 'MarkerSize', 8+abl*3, ...
                'MarkerFaceColor', col(ild,:), 'MarkerEdgeColor', 'w')
        end
    end
    
    mu = mean(mean(nanmean(XX,4),3),2);
    plot3([0 mu(1)], [0 mu(2)], [0 mu(3)], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    set(ann, 'string', labels{ii})
    
    thetas = 1:2:360;
    for frame = 1:length(thetas)
        view([cosd(-135+thetas(frame)) sind(-135+thetas(frame)) 0.2])
        pause(0.1)
        framesize = 44;
        thisFrame = getframe(gcf,[(1000-framesize*16)/2 (1000-framesize*9)/2 framesize*16 framesize*9]);
        writeVideo(vidObj, thisFrame);
    end
end

close(vidObj)

    function U = gramschmidt(V)
        n = size(V,1);
        k = size(V,2);
        U = zeros(n,k);
        U(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
        for i = 2:k
            U(:,i) = V(:,i);
            for j = 1:i-1
                U(:,i) = U(:,i) - ( U(:,i)'*U(:,j) )/( U(:,j)'*U(:,j) )*U(:,j);
            end
            U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
        end
    end


    function h = plotEllipsoid(data)
        scale = sqrt(chi2inv(0.1, 3));

        data = data(~isnan(sum(data,2)), 1:3);
        data0 = bsxfun(@minus, data, mean(data));

        [v, d] = eig(cov(data0));
        [d, ind] = sort(diag(d), 'descend');
        d = diag(d);
        v = v(:, ind);

        t = linspace(0,2*pi,32);
        theta = linspace(-pi,pi,32)';

        e = [];
        for i=1:length(theta)
            e = [e [cos(t)*cos(theta(i)); sin(t)*cos(theta(i)); ones(size(t))*sin(theta(i))]];
        end
        vv = v*sqrt(d) * scale;               
        e = bsxfun(@plus, vv*e, mean(data)');
        h = surf(reshape(e(1,:),32,[]), reshape(e(2,:),32,[]), reshape(e(3,:),32,[]));
    end
end
