%% Calcium Imaging
%% (1) self (pre) configuration 
DATA_DIR = pwd;
% Change the directory to where the data set are 

preExp = struct();
preExp.Names = {dir(fullfile(DATA_DIR,'*-pre.mat')).name};
preExp.Names = cellfun(@(x) x(1:7),preExp.Names,'UniformOutput',false); %for WT (mPFCXXX)
%preExp.Names = cellfun(@(x) x(1:13),preExp.Names,'UniformOutput',false); %for Shank3-KO (shank3mPFCXXX)
%preExp.Names = cellfun(@(x) x(1:6),preExp.Names,'UniformOutput',false); %for Shank3-cKO (cKOXXX)
preExp.nRecordings = length(preExp.Names);

PreShock = struct();
PreShock.Delay = 140;
PreShock.Duration = 2;
PreShock.Interval = 0;
PreShock.nShocks = 1;
PreShock.Onsets = (0:PreShock.nShocks-1)*(PreShock.Duration+PreShock.Interval) + PreShock.Delay;

%% Import raw PreShockCa transient data, calculate zscored df/f 
PreShockCa = struct();
PreShockCa.Names = {dir(fullfile(DATA_DIR,'*-pre.mat')).name};
PreShockCa.Names = cellfun(@(x) x(1:7),preExp.Names,'UniformOutput',false); % For mPFC (mPFC00X)
%PreShockCa.Names = cellfun(@(x) x(1:13),preExp.Names,'UniformOutput',false); % For Shank3-KO (shank3mPFC00X)
%PreShockCa.Names = cellfun(@(x) x(1:6),preExp.Names,'UniformOutput',false); % For Shank3-cKO (cKO00X)

for ii = 1:preExp.nRecordings
    ct = [];
    load([preExp.Names{ii} '-pre.mat'],'output');
    T = output.temporal_weights;
    labels = output.labels;
    numrow = length(labels(1,:));
    for k = 1:1:numrow 
        if labels(:,k) == 1
           ct = [ct (T(:,k))];
        end
    end 
    PreShockCa.dff{ii} = ct';
    PreShockCa.dffzs{ii} = zscore(ct',0,2);
    clear T labels numrow a
end

clear numel output ct zs dfm ii k 

%% Get df/f 
PreShockCa.shock = struct();
PreShockCa.shock.t = 2; % 
sec = PreShockCa.shock.t;
for ii = 1:preExp.nRecordings
    PreShockCa.shock.BA{ii} = PreShockCa.dffzs{ii}(:,(140-sec)*25+1:(140+sec)*25); % -10 ~ +10
end
clear ii j 
clear k 

%% Permutation test for identifying the self-shock responding neurons (cirshift)
sec = PreShockCa.shock.t; 
for ii = 1:preExp.nRecordings
    num_cells = size(PreShockCa.dffzs{ii}, 1);
    num_samples_per_session = 25 * 180; % 25 Hz for 180 seconds

    shockA = PreShockCa.dffzs{ii}(:, (140-sec)*25+1:140*25);
    shockB = PreShockCa.dffzs{ii}(:, 140*25:(140+sec)*25);

    % Num of shuf iterations
    num_iterations = 10000;

    % Calculate the actual difference between means
    actual_diff = mean(shockB, 2) - mean(shockA, 2);
    PreShockCa.shock.actual_diff{ii} = actual_diff;
    [~, PreShockCa.shock.order{ii}] = sort(PreShockCa.shock.actual_diff{ii});

    % Perform circular shifting between shock A and shock B
    for i = 1:num_cells
        shifted_diff = zeros(1, num_iterations);
        for j = 1:num_iterations
            % Circularly shift the data between the two sessions
            shift_amount = randi([1, num_samples_per_session]); % random shift amount
            shifted_preshock = circshift(PreShockCa.dffzs{ii}(i, :), shift_amount,2);% Circular shift along the second dimension
            shifted_A = shifted_preshock(:,1:sec*25);
            shifted_B = shifted_preshock(:,sec*25+1:sec*2*25);
            % Calculate the difference between means for shifted data
            shifted_diff(j) = mean(shifted_B) - mean(shifted_A);
        end

        % Calculate the 2.5% and 97.5% quantiles of the shifted differences
        quantile_2_5 = quantile(shifted_diff, 0.025);
        quantile_97_5 = quantile(shifted_diff, 0.975);

        % Check if the actual difference falls within the quantiles
        if actual_diff(i) > quantile_97_5
            PreShockCa.shock.IDcir{ii}(i) = 1;
        elseif actual_diff(i) < quantile_2_5
            PreShockCa.shock.IDcir{ii}(i) = -1;
        else
            PreShockCa.shock.IDcir{ii}(i) = 0;
        end
    end
end


clear i ii num_* shifted_* j actual_diff
clear deltaFF_diff quantile_* sec shock* k
clear shift_amount

%% Fig.2fgh
a = []; b = []; c = [];
for ii = 1:preExp.nRecordings
    Ca = PreShockCa.dffzs;
    labels = PreShockCa.shock.IDcir;
    w = size(PreShockCa.dffzs{ii},1);
    for k = 1:w
        if labels{ii}(k) == 1
           a = [a; (Ca{ii}(k,3251:3750))];
        elseif labels{ii}(k) == -1
                b = [b; (Ca{ii}(k,3251:3750))];
        elseif labels{ii}(k) == 0
                c = [c; (Ca{ii}(k,3251:3750))];
        end 
    end
end

figure 
x = 1:size(a,2);
hold on 
shadedErrorBar(x,mean(a,1,'omitnan'),std(a,'omitnan')/sqrt(size(a,1)),'lineProps','g');
shadedErrorBar(x,mean(b,1,'omitnan'),std(b,'omitnan')/sqrt(size(b,1)),'lineProps','y');
shadedErrorBar(x,mean(c,1,'omitnan'),std(c,'omitnan')/sqrt(size(c,1)),'lineProps','k');
ylim([-1.0 5.5]) 
hold off

aa = mean(a,'omitnan')'; asem = std(a,'omitnan')'/sqrt(length(a))';
bb = mean(b,'omitnan')'; bsem = std(b,'omitnan')'/sqrt(length(b))';
cc = mean(c,'omitnan')'; csem = std(c,'omitnan')'/sqrt(length(c))';
PreShockCa.ShockRescir.trace = [aa,asem,bb,bsem,cc,csem];

clear aa bb cc dd ee ff asem bsem csem dsem esem fsem k a b c Ca d e f ii v w x

%% (2) OF import & analysis 
DATA_DIR = pwd;
% Change the directory to where the data set are 


pallet = get(groot,'DefaultAxesColorOrder');

Exp = struct();
Exp.Names = {dir(fullfile(DATA_DIR,'*_EXTRACT.mat')).name};
Exp.Names = cellfun(@(x) x(1:7),Exp.Names,'UniformOutput',false); %for WT (mPFCXXX)
%Exp.Names = cellfun(@(x) x(1:13),Exp.Names,'UniformOutput',false); %for Shank3-KO (shank3mPFCXXX)
%Exp.Names = cellfun(@(x) x(1:6),Exp.Names,'UniformOutput',false); %for Shank3-cKO (cKOXXX)

% 
Exp.nRecordings = length(Exp.Names);
Exp.TrialLength = 900; % sec
Exp.nCells = zeros(1,Exp.nRecordings);

Shock = struct();
Shock.Delay = 300;
Shock.Duration = 2;
Shock.Interval = 8;
Shock.nShocks = 60;
Shock.Onsets = (0:Shock.nShocks-1)*(Shock.Duration+Shock.Interval) + Shock.Delay;


%% Load calcium event data with OASIS

Event = struct();
Event.Fs = 25;

Event.FrameLength = 1/Event.Fs;
Event.ChunkLength = 2;
Event.EpochLength = 10; % sec

Event.nFrames = Exp.TrialLength / Event.FrameLength;
Event.nChunks = Exp.TrialLength / Event.ChunkLength;
Event.nEpochs = Exp.TrialLength / Event.EpochLength;

Event.nFramesPerChunk = Event.ChunkLength / Event.FrameLength;
Event.nChunksPerEpoch = Event.EpochLength / Event.ChunkLength;
Event.nFramesPerEpoch = Event.EpochLength / Event.FrameLength;

Event.tFrame = 0:Event.FrameLength:Exp.TrialLength;
Event.tChunk = 0:Event.ChunkLength:Exp.TrialLength;
Event.tEpoch = 0:Event.EpochLength:Exp.TrialLength;

Event.Data = cell(Exp.nRecordings,1);
Event.eventTimes = cell(Exp.nRecordings,1);
Event.nEvents = struct();
        
for ii=1:Exp.nRecordings
    Event.Data{ii} = ...
        load(fullfile(DATA_DIR,[Exp.Names{ii} '_OASIS.mat'])).S; % or S_tmp
    Event.Data{ii}(Event.Data{ii}>0) = 1;
    Exp.nCells(ii) = size(Event.Data{ii,:},2);
    for jj=1:Exp.nCells(ii)
            myn = Event.Data{ii}(:,jj);
            if size(myn,1) < 22500
                myn = [myn;NaN(22500-size(myn,1),1)];
            elseif size(myn,1) == 22500
                myn = myn;
            elseif size(myn,1) > 22500
                myn = myn(1:22500,:);
            end
        for period = ["Chunk","Epoch"]
            d1 = Event.(strcat("nFramesPer",period));
            d2 = Event.(strcat("n",period,"s"));
            Event.nEvents(ii).(strcat("By",period))(:,jj) = sum(reshape(myn,d1,d2));
        end
    end
end

clear ii jj my* d1 d2 period

%% Import raw calcium transient data  

Calcium = struct();
for ii = 1:Exp.nRecordings
    ct = [];
    load([Exp.Names{ii} '_EXTRACT.mat'],'output');
    %load([Exp.Names{ii} '_EXTRACT_output.mat'],'output');
    T = output.temporal_weights;
    labels = output.labels;
    numrow = length(labels(1,:));

    for k = 1:1:numrow 
        if labels(:,k) == 1
           ct = [ct (T(:,k))];
        end
    end
 
    Calcium.dff{ii} = ct';
    Calcium.dffzs{ii} = zscore(ct',0,2);
    if size(Calcium.dffzs{ii},2) == 22500
        Calcium.dffzs{ii} = Calcium.dffzs{ii};
    elseif size(Calcium.dffzs{ii},2) < 22500
        a = 22500 - size(Calcium.dffzs{ii},2);
        Calcium.dffzs{ii} = [Calcium.dffzs{ii}, NaN(size(Calcium.dffzs{ii},1),a)];
    elseif size(Calcium.dffzs{ii},2) > 22500
        Calcium.dffzs{ii} = Calcium.dffzs{ii}(:,1:22500);
    end
    clear T labels numrow a
end

clear numel output ct zs dfm ii k 

%% Alignment and average the ca2+ trace (conditioning 301s~900s) 

for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
        for j = 1:1:w
            E = [];
            for k = 7451:250:22301
                E = [E; Calcium.dffzs{ii}(j,k:k+299)];
            end
            Calcium.RawShockDffzs{ii}{j} = E;
        end
end

for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    for j = 1:1:w
       Calcium.RawShockDffzsMean{ii}{j} = mean(Calcium.RawShockDffzs{ii}{j},'omitnan');
    end
end

clear ii w j D E k

% Alignment and average the ca2+ trace (habituation 0~300s) 

for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    
    for j = 1:1:w
        E = [];
        for k = 201:250:7201
            E = [E; Calcium.dffzs{ii}(j,k:k+299)];
        end
        Calcium.RawShockDffzsHabi{ii}{j} = E;
    end
end

for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    for j = 1:1:w
        Calcium.RawShockDffzsMeanHabi{ii}{j} = mean(Calcium.RawShockDffzsHabi{ii}{j},'omitnan');
    end
end
clear ii w j D E k w

%% Calculate the mean of the trace (-2~0, 0~2) for shock responding neurons
% (Shock responding neurons' annotation) 
for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    for j = 1:w
        A = [];
        B = [];       
        for k = 7451:250:22250  % Shock (conditioning) period 
            A = [A; mean(Calcium.dffzs{ii}(j,k:k+49))];
            B = [B; mean(Calcium.dffzs{ii}(j,k+50:k+99))];
            % B = [B; mean(Calcium.dffzs{ii}(j,k+50:k+99)),omitnan];
        end
        a = [A, B];
        Calcium.AvgShockdffzs{ii}{j} = a;
        clear a
    end
end
clear A B C D

% mean of the trace
for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    for j = 1:w
        Calcium.ConddffzsMean{ii}{j} = mean(Calcium.AvgShockdffzs{ii}{j},'omitnan');
    end
end
clear w v j ii k

% Stats
Calcium.ttest.dffzs = struct();
for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    pabc = NaN(w,1);
    habc = NaN(w,1);  
    for j = 1:w
        a = Calcium.AvgShockdffzs{ii}{j}(:,1);
        b = Calcium.AvgShockdffzs{ii}{j}(:,2);      
        [ha,pa] = ttest(a,b); %for comparing -2s~0s vs. 0s~2s      
        pabc(j,:) = pa; 
        habc(j,:) = ha;
    end
    Calcium.ttest.dffzs.pvalue{ii} = pabc;
    Calcium.ttest.dffzs.h{ii} = habc;
    
end
clear pabc habc ii w v a b j k ha pa 

% sign with df/f
for ii = 1:Exp.nRecordings
    w = size(Calcium.dffzs{ii},1);
    for j = 1:w
        Calcium.ShocksignCond.dffzs{ii}(j,1) = sign(Calcium.ConddffzsMean{ii}{j}(1,2) - Calcium.ConddffzsMean{ii}{j}(1,1));
    end 
end
clear w ii j k 

% annotation for 0-2 vs. -2-0 
Calcium.annotation.dffzs = struct();
Calcium.annotation.dffzs.GroupNames = {'NegMod', 'PosMod'};
Calcium.annotation.dffzs.GroupFlags = [-1,1];
Calcium.annotation.dffzs.Thresh = 0.05;

Bb = []; Cc = []; Dd = [];
for ii = 1:Exp.nRecordings
    Calcium.annotation.dffzs.group{ii} = Calcium.ShocksignCond.dffzs{ii} .* (Calcium.ttest.dffzs.pvalue{ii} < Calcium.annotation.dffzs.Thresh);
    Bb = [Bb; sum(Calcium.annotation.dffzs.group{ii}>0)];
    Cc = [Cc; sum(Calcium.annotation.dffzs.group{ii}<0)];
    Dd = [Dd; sum(Calcium.annotation.dffzs.group{ii}==0)];
end


Calcium.annotation.dffzs.NumPosMod = Bb;
Calcium.annotation.dffzs.NumNegMod = Cc;
Calcium.annotation.dffzs.Numns = Dd; 

% # of Shock act; n.s.; Shock sup.
Calcium.annotation.dffzs.NumAll = [Bb(:,1), Dd(:,1), Cc(:,1)];
Calcium.annotation.dffzs.NumAllSum = sum(Calcium.annotation.dffzs.NumAll);
Calcium.annotation.dffzs.NumAllSum

clear Bb Cc Dd Cabc ii ans

%% Fig.2ijk
a = []; b = []; c = [];

for ii = 1:Exp.nRecordings
    Ca = Calcium.RawShockDffzs;
    labels = Calcium.annotation.dffzs.group;
    w = size(Calcium.dffzs{ii},1);
    for k = 1:w
        if labels{ii}(k) == 1
           a = [a; (Ca{ii}{k})];
        else if labels{ii}(k) == -1
                b = [b; (Ca{ii}{k})];
        else 
                c = [c; (Ca{ii}{k})];
            end 
        end
    end
end

clear Ca 


d = []; e = []; f = [];
for ii = 1:Exp.nRecordings
    Ca = Calcium.RawShockDffzsHabi;
    labels = Calcium.annotation.dffzs.group;
    w = size(Calcium.dffzs{ii},1);

    for k = 1:w
        if labels{ii}(k) == 1
           d = [d; (Ca{ii}{k})];
        else if labels{ii}(k) == -1
                e = [e; (Ca{ii}{k})];
        else 
                f = [f; (Ca{ii}{k})];
            end 
        end
    end
end

%
figure %habituation
x = 1:size(a,2);
hold on 
shadedErrorBar(x,mean(d,1,'omitnan'),std(d,'omitnan')/sqrt(size(d,1)),'lineProps','g');
shadedErrorBar(x,mean(e,1,'omitnan'),std(e,'omitnan')/sqrt(size(e,1)),'lineProps','y');
shadedErrorBar(x,mean(f,1,'omitnan'),std(f,'omitnan')/sqrt(size(f,1)),'lineProps','k');
ylim([-0.3 0.3]) 
hold off

figure %conditioning
x = 1:size(a,2);
hold on 
shadedErrorBar(x,mean(a,1,'omitnan'),std(a,'omitnan')/sqrt(size(a,1)),'lineProps','g');
shadedErrorBar(x,mean(b,1,'omitnan'),std(b,'omitnan')/sqrt(size(b,1)),'lineProps','y');
shadedErrorBar(x,mean(c,1,'omitnan'),std(c,'omitnan')/sqrt(size(c,1)),'lineProps','k');
ylim([-0.3 0.3])
hold off

aa = mean(a,'omitnan')'; asem = std(a,'omitnan')'/sqrt(length(a))';
bb = mean(b,'omitnan')'; bsem = std(b,'omitnan')'/sqrt(length(b))';
cc = mean(c,'omitnan')'; csem = std(c,'omitnan')'/sqrt(length(c))';
dd = mean(d,'omitnan')'; dsem = std(d,'omitnan')'/sqrt(length(d))';
ee = mean(e,'omitnan')'; esem = std(e,'omitnan')'/sqrt(length(e))';
ff = mean(f,'omitnan')'; fsem = std(f,'omitnan')'/sqrt(length(f))';
Calcium.ShockRes.trace = [dd,dsem,ee,esem,ff,fsem,aa,asem,bb,bsem,cc,csem];

clear aa bb cc dd ee ff asem bsem csem dsem esem fsem k a b c Ca d e f ii v w x


%% Shock Event cosine similarity 

for ii=1:Exp.nRecordings
    myf = repmat([1 0 0 0 0],1,60)'; % 2s shock 
    myn = Event.nEvents(ii).ByChunk(151:450,:);
    for q = 1:size(Event.nEvents(ii).ByChunk(151:450,:),2)
        A = [myf,myn(:,q)];
        A(any(isnan(A),2),:)= []; % handling NaN value 
        Event.ShockPeriodCosSim.BCr{ii}(q) = getCosineSimilarity(A(:,1),A(:,2));
        clear A 
    end 
end

clear ii jj my* q clear A 

%% (3) Match (pre + OF) 
%% Collect cellregID
CellRegID = struct();
CellRegID.Names = Exp.Names; 

for ii = 1:Exp.nRecordings
    load([Exp.Names{ii} '_CellReg.mat'],'cell_registered_struct');
    Tc = cell_registered_struct.cell_to_index_map;
    T = Tc; % becasue 1st column is OF, 2nd column is pre
    numZero = sum(T(:,1) == 0);
    Ta = sortrows(T);
    Tb = [Ta(numZero+1:end,:);Ta(1:numZero,:)];
    CellRegID.Index{ii} = Tb; % make OF-pre order 
    clear T Ta Tb numZero Tc
end
clear ii cell_registered_struct

%% Collect 
CellRegID.PreShockMatch = [];
for ii =1:Exp.nRecordings
        CellRegID.PreShockMatchcir{ii} = CellRegID.Index{ii};
    for j = 1:length(CellRegID.Index{ii})
        if CellRegID.Index{ii}(j,2) == 0
            CellRegID.PreShockMatchcir{ii}(j,3) = NaN;
        elseif CellRegID.Index{ii}(j,2) ~= 0 && PreShockCa.shock.IDcir{ii}(CellRegID.Index{ii}(j,2)) == 1
            CellRegID.PreShockMatchcir{ii}(j,3) = 1;
        elseif CellRegID.Index{ii}(j,2) ~= 0 && PreShockCa.shock.IDcir{ii}(CellRegID.Index{ii}(j,2)) == -1
            CellRegID.PreShockMatchcir{ii}(j,3) = -1;
        else
            CellRegID.PreShockMatchcir{ii}(j,3) = 0; 
        end   
    end
     
end
clear ii j k 

%% Labels 

for k = 1:3
    for ii = 1:Exp.nRecordings
        a = sum(CellRegID.PreShockMatchcir{ii}(:,1)~=0);
        CellRegID.labelscir{ii} = CellRegID.PreShockMatchcir{ii}(1:a,3);
    end
end

clear a ii k 

%% Alignment 
%#1:preshock num, %2:OF num, #3: ShockCosSim (OF),
%#4: pre-actual diff, #5: Pre-shock id(circshift) 

for ii = 1:Exp.nRecordings
    b = Event.ShockPeriodCosSim.BCr{ii}';
    s = size(CellRegID.Index{ii},1)-size(b,1);
    CellRegID.align.DifsCcS{ii} = [CellRegID.Index{ii}(:,2),CellRegID.Index{ii}(:,1),[b;NaN(s,1)]];
    % align with the OF-tracked cell first, coossim
    CellRegID.align.DifsCcS{ii} = sortrows(CellRegID.align.DifsCcS{ii},1);
    % align with pre-shock cells 
    sn = sum(CellRegID.align.DifsCcS{ii}(:,1) == 0);%zeros of the preshock neurons
    % delete the un-tracked cells 
    c = [CellRegID.align.DifsCcS{ii}(sn+1:end,:),PreShockCa.shock.actual_diff{ii},PreShockCa.shock.IDcir{ii}'];
   
    CellRegID.align.DifsCcS{ii} = c;
end

clear q ii a b d s sn c k 

% Collect for all the cells 
X = [];
for ii = 1:Exp.nRecordings
    X = [X; CellRegID.align.DifsCcS{ii}];
end
CellRegID.align.summary = X;
clear X q ii 


%% Fig.2rst visualization

Z = CellRegID.align.summary;
data = [Z(:,3),Z(:,4)]; % OF-Shock cossim vs. pre-shock act dif

% for graph
x = data(:, 1);
y = data(:, 2);
[myr, myp] = corr(x, y, 'Rows', 'pairwise');

% Remove rows with NaN
validIdx = ~isnan(x) & ~isnan(y);
x = x(validIdx);
y = y(validIdx);

% Fit a linear regression model
X = [ones(size(x)) x];
b = X\y;
yfit = X * b;
slope = b(2); 

% Plot
figure;
scatter(x, y, 'filled');
hold on;
plot(x, yfit, 'r-', 'LineWidth', 1);
title(sprintf('Correlation coefficient: %.2f (p = %.2g)', myr, myp));
legend('Data', 'Fit', 'Location', 'Best');
xlim([0 0.7])
ylim([-4 8])
pbaspect([1 1 1]); 
hold off;

clear x y my* b yfit X Z 

%% Functions
function Cs = getCosineSimilarity(x,y)
% 
% call:
% 
%      Cs = getCosineSimilarity(x,y)
%      
% Compute Cosine Similarity between vectors x and y.
% x and y have to be of same length. The interpretation of 
% cosine similarity is analogous to that of a Pearson Correlation
% 
% R.G. Bettinardi
% -----------------------------------------------------------------


if isvector(x)==0 || isvector(y)==0
    error('x and y have to be vectors!')
end

if length(x)~=length(y)
    error('x and y have to be same length!')
end

xy   = dot(x,y);
nx   = norm(x);
ny   = norm(y);
nxny = nx*ny;
Cs   = xy/nxny;
end
