%% Load data 

load('240521-nWT-OASIS-3.mat')
% load('240524-KO_OASIS.mat')

Event.PreData{1} = load('oasis_WT001_pre').S_tmp;
Event.PreData{2} = load('oasis_WT002_pre').S_tmp;
Event.PreData{3} = load('oasis_WT003_pre').S_tmp;
Event.PreData{4} = load('oasis_WT004_pre').S_tmp;
% Event.PreData{1} = load('oasis_KO001_pre').S_tmp;
% Event.PreData{2} = load('oasis_KO006_pre').S_tmp;
% Event.PreData{3} = load('oasis_KO008_pre').S_tmp;
% Event.PreData{4} = load('oasis_KO009_pre').S_tmp;
Event.PreData{1}(Event.PreData{1}>0) = 1;
Event.PreData{2}(Event.PreData{2}>0) = 1;
Event.PreData{3}(Event.PreData{3}>0) = 1;
Event.PreData{4}(Event.PreData{4}>0) = 1;

%% Individual PCA + Mahalanobis distances

id1 = find(CellRegID.align.summary{1,1}(:,1) ==1);
% CellRegID.align.summary(:,1) ~=0 && CellRegID.align.summary(:,2) ~=0

% meanAllSO = zeros(1000,4); % self-other
meanAllON = zeros(1000,4); % other-no
% meanAllSS = zeros(1,4); % self-self
meanAllOO = zeros(1,4); % other-other

% num of dimensions to use to calculate mahal. dis.
nDim = 7;

selfShock = 701:710; % 1:900 self-shock period
otherShock = [];
for startVal = 2401:50:5351 % 2401:5400 Cond(OF) period
    otherShock = [otherShock startVal:startVal+9];
end
noShock = 1:700; 

for nn=1:4
    X = [];
    % X = Event.nEvents(nn).ByChunk';
    for n = 1:size(cell2mat(Event.PreData(nn)),2)
        m = n +id1(nn)-1;
        if CellRegID.align.summary{1,1}(m,2) >0
            line = [Event.PreData{nn}(:,n); Event.Data{nn}(:,CellRegID.align.summary{1,1}(m,2))];
            % line = Event.Data{nn}(:,CellRegID.align.summary{1,1}(m,2));
            X = [X; line'];
        end
    end
    
    % 5Hz
    X_2 = zeros(size(X,1), round(size(X,2)/5));
    for t = 1:size(X,2)/5
        for n =1:size(X,1)
            X_2(n,t) = sum(X(n,t*5-4:t*5));
        end
    end
    
    [coeff,score,latent,tsquared,explained] = pca(X_2(:,1:end));
    
    % Fig.3a
    figure
    plot(coeff(1:700,1),coeff(1:700,2),'k','LineWidth',0.1)
    hold on
    for n = 1:60
        scatter(coeff(2401+(n-1)*50:2410+(n-1)*50,1),coeff(2401+(n-1)*50:2410+(n-1)*50,2),15,'b','filled')
        plot(coeff(2401+(n-1)*50:2410+(n-1)*50,1),coeff(2401+(n-1)*50:2410+(n-1)*50,2),'b','LineWidth', 0.2)
        hold on
    end
    scatter(coeff(701:710,1),coeff(701:710,2),15,'r', 'filled')
    plot(coeff(701:710,1),coeff(701:710,2),'r','LineWidth',0.2)
    hold on
    ax = gca;
    ax.YLim = [-0.04 0.07];
    ax.PlotBoxAspectRatio = [1 1 1];

    % Calculate Mahalanobis Distance
    dataSelfShock = coeff(selfShock,1:nDim); % 10
    dataOtherShock = coeff(otherShock,1:nDim); % 600
    dataNoShock = coeff(noShock,1:nDim); % 1500
    
    covMatrix = cov(coeff(:,1:nDim));
    
    % meanSO = [];
    meanON = [];
    rng default  
    % data for Fig.3c
    for r = 1:1000
        tmp = randi([1,600],10,1);
        dataA = dataOtherShock(tmp,:);
        % dataB = dataSelfShock;
        tmp2 = randi([1,700],10,1);
        dataB = dataNoShock(tmp2,:);
    
        distanceAll = [];
        for n = 1:size(dataA,1)
            x = dataA(n, :);
            for m = 1:size(dataB,1)
                y = dataB(m, :);
                distance = mahalanobisDistance(x, y, covMatrix);
                distanceAll = [distanceAll; distance];
            end
        end
        % meanSO = [meanSO; mean(distanceAll)];
        meanON = [meanON; mean(distanceAll)];
    end
    % meanAllSO(:,nn) = meanSO;
    meanAllON(:,nn) = meanON;
    
    % meanSS = [];
    meanOO = [];
    % data for Fig.3d
    % dataA = dataSelfShock;
    dataA = dataOtherShock;

    distanceAll = [];
    for n = 1:size(dataA,1)
        x = dataA(n, :);
        for m = n+1:size(dataA,1)
            y = dataA(m, :);
            distance = mahalanobisDistance(x, y, covMatrix);
            distanceAll = [distanceAll; distance];
        end
    end
    % meanSS = [meanSS; mean(distanceAll)];
    meanOO = [meanOO; mean(distanceAll)];
    % meanAllSS(:,nn) = meanSS;
    meanAllOO(:,nn) = meanOO;
end

%% Function
function d = mahalanobisDistance(x, y, covMatrix)
    delta = x - y;
    d = sqrt(delta * (covMatrix \ delta'));
end