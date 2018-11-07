function makeMovies(l)

% Load the data if not provided
if nargin==0
    l = load('evokedResponses_150to150.mat');
end

begtime = 0.00;
endtime = 0.15;

example_ind = [1 7];
col = parula12();
filenames = {'movie_inactive.avi', 'movie_active.avi'};

for ii = 1
    f = figure('Position', [0 0 800 800]); 
    set(gca,'Color','k')    
    axis(110*[-1 1 -1 1 -1 1])
    set(gcf,'Renderer','zbuffer') 
    axis normal
    hold on
    axis vis3d
    
    d = example_ind(ii);    
    X = squeeze(sum(l.datasets{d}(:,:,:,l.time>begtime & l.time<endtime ,:),4) / (endtime-begtime));
    XX = X(:,:)';
    
    Xpsth = nanmean(X,4);
    mu = mean(Xpsth(:,:),2);
    Xabl = squeeze(mean(Xpsth,2))';
    [~,~,Vabl] = svd(Xabl - mean(Xabl), 'econ');
    Xild = squeeze(nanmean(Xpsth,3))';
    [~,~,Vild] = svd(Xild - mean(Xild), 'econ');
    V = gramschmidt([mu Vild(:,1) Vabl(:,1)]);
    
    XX  = reshape((XX * V)', 3, 12, 3, []);

    % ellipsoids
    for ild=1:12
        for abl=1:3
            data = squeeze(XX(:,ild,abl,:))';
            h = plotEllipsoid(data);
            set(h, 'LineStyle', 'none')
            set(h, 'FaceColor', col(ild,:))
            set(h, 'FaceAlpha', 0.03)
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
    
    % zero
    plot3(0, 0, 0, 'o', 'MarkerSize', 8+2*3, ...
        'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'w')
    mu = mean(mean(nanmean(XX,4),3),2);
    plot3([0 mu(1)], [0 mu(2)], [0 mu(3)], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    % make movie
%     h = rotate3d(gca());
%     set(h, 'Enable', 'on')

    vidObj = VideoWriter('test.avi');
    open(vidObj);
    thetas = 1:10:360;
    for frame = 1:length(thetas)
        view([cosd(-135+thetas(frame)) sind(-135+thetas(frame)) 0.2])
        pause(0.01)
        thisFrame = getframe(gcf,[0 0 700 700]);
        thisFrame.cdata = thisFrame.cdata(1:700,1:700,:);
        writeVideo(vidObj, thisFrame);
    end
    close(vidObj);
end

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
