%% ========================================================================
%  SFig4_analysis.m
%
%  Supplementary Figure 4
%  Shock-locked center-position dynamics of behavioural components
%  in observer mice.
%
%  Input arrays are outputs of BehavClassification.m, which tracks observer
%  mice with DeepLabCut and assigns each 2-s bout to a behavioural cluster
%  by t-SNE embedding followed by watershed segmentation.
%
%  Panels
%  ------
%  (a) Mean shock-aligned displacement (dX) over time, per trial group.
%  (b) Mean dX averaged over 0-10 s after shock onset, per trial group.
%  (c) Signed area under the dX curve over the same window, per trial group.
%
% ========================================================================

clear; close all; clc;

load('SFig4_data.mat');   % cluWTRaw2s, bout_x

%% ------------------------------------------------------------------------
%  Parameters
%  ------------------------------------------------------------------------
% Behavioural components (raw cluster indices)
i1i2   = [20 18];        % turned-away / escape component
i3i4i5 = [23 27 25];     % other escape-like immobile components

% Session structure
nAnimals    = size(cluWTRaw2s, 1);
nBouts      = size(cluWTRaw2s, 2);
boutDur     = 2;                    % s per bout
shockBouts  = 151:5:446;            % 60 foot shocks, one every 10 s
nShocks     = numel(shockBouts);

% Peri-shock window, in bouts before and after the shock bout
preBouts    = 4;                    % window starts at -8 s
postBouts   = 4;                    % window ends   at +8 s

% Analysis windows, defined on the bout start time (s from shock onset)
baseWin     = [-8 -2];              % baseline: bouts spanning -8 to 0 s
classWin    = [-2  6];              % component must occur within this window
postWin     = [ 0  8];              % bouts spanning 0 to 10 s, used for b and c

% Plotting
markerShift = 1;                    % plot each bout at its centre (+1 s)
colGroup    = {[0.90 0.49 0.13], ...   % i1+i2      orange
               [0.58 0.40 0.74], ...   % i3+i4+i5   purple
               [0.40 0.40 0.40]};      % rest       grey

%% ------------------------------------------------------------------------
%  Build the peri-shock window
%  ------------------------------------------------------------------------
nWin    = preBouts + postBouts + 1;
timeAx  = (-preBouts:postBouts) * boutDur;      % bout start time
timePlot = timeAx + markerShift;                % x-coordinate for plotting

baseBins  = find(timeAx >= baseWin(1)  & timeAx <= baseWin(2));
classBins = find(timeAx >= classWin(1) & timeAx <= classWin(2));
postBins  = find(timeAx >= postWin(1)  & timeAx <= postWin(2));

fprintf('Peri-shock window : %s s\n', mat2str(timeAx));
fprintf('Baseline bouts    : %s s\n', mat2str(timeAx(baseBins)));
fprintf('Grouping bouts    : %s s\n', mat2str(timeAx(classBins)));
fprintf('Post-shock bouts  : %s s\n\n', mat2str(timeAx(postBins)));

%% ------------------------------------------------------------------------
%  Collect x position and cluster identity for every trial
%  ------------------------------------------------------------------------
xPeri   = nan(nAnimals, nShocks, nWin);
cluPeri = zeros(nAnimals, nShocks, nWin);

for a = 1:nAnimals
    for s = 1:nShocks
        for w = 1:nWin
            b = shockBouts(s) - preBouts + w - 1;
            if b < 1 || b > nBouts, continue; end
            xPeri(a,s,w)   = bout_x(a,b);
            cluPeri(a,s,w) = cluWTRaw2s(a,b);
        end
    end
end

% Subtract each trial's own pre-shock baseline
baseline = nanmean(xPeri(:,:,baseBins), 3);
dxPeri   = xPeri - baseline;

% Flatten to one row per trial
dx = reshape(dxPeri, nAnimals*nShocks, nWin);

%% ------------------------------------------------------------------------
%  Group trials by behavioural component
cluClass = cluPeri(:,:,classBins);

hasI1I2   = squeeze(any(ismember(cluClass, i1i2),   3));
hasI3I4I5 = squeeze(any(ismember(cluClass, i3i4i5), 3));

isI1I2   = hasI1I2(:);
isI3I4I5 = hasI3I4I5(:);
isRest   = ~isI1I2 & ~isI3I4I5;

trialGroup = {isI1I2, isI3I4I5, isRest};
groupName  = {'i1+i2', 'i3+i4+i5', 'rest'};
nGroup     = 3;

fprintf('Trial counts\n');
for g = 1:nGroup
    fprintf('  %-10s n = %4d\n', groupName{g}, sum(trialGroup{g}));
end
fprintf('  %-10s n = %4d\n', 'both', sum(isI1I2 & isI3I4I5));
fprintf('  %-10s n = %4d\n\n', 'total', nAnimals*nShocks);

%% ------------------------------------------------------------------------
%  Panel a: mean +/- SEM time course
%  ------------------------------------------------------------------------
dxMean = nan(nGroup, nWin);
dxSEM  = nan(nGroup, nWin);
groupN = nan(1, nGroup);

for g = 1:nGroup
    d = dx(trialGroup{g}, :);
    groupN(g)   = size(d, 1);
    dxMean(g,:) = nanmean(d, 1);
    dxSEM(g,:)  = nanstd(d, 0, 1) ./ sqrt(sum(~isnan(d), 1));
end

%% ------------------------------------------------------------------------
%  Panels b and c: per-trial summary measures over 0-10 s
%  ------------------------------------------------------------------------
tPost   = timeAx(postBins);
meanDX  = cell(1, nGroup);
aucDX   = cell(1, nGroup);

for g = 1:nGroup
    d = dx(trialGroup{g}, postBins);
    meanDX{g} = nanmean(d, 2);

    auc = nan(size(d,1), 1);
    for t = 1:size(d,1)
        v  = d(t,:);
        ok = ~isnan(v);
        if nnz(ok) >= 2
            auc(t) = trapz(tPost(ok), v(ok));
        end
    end
    aucDX{g} = auc;
end

%% ------------------------------------------------------------------------
%  Statistics
%  ------------------------------------------------------------------------
pMeanDX = nan(1, nGroup);
pAUC    = nan(1, nGroup);

fprintf('One-sample two-sided t-test against baseline\n');
fprintf('%-10s %6s %10s %8s %12s %10s %10s %12s\n', ...
        'Group', 'n', 'mean dX', 'SEM', 'p (dX)', 'mean AUC', 'SEM', 'p (AUC)');

for g = 1:nGroup
    b = meanDX{g}(~isnan(meanDX{g}));
    c = aucDX{g}(~isnan(aucDX{g}));

    [~, pMeanDX(g)] = ttest(b);
    [~, pAUC(g)]    = ttest(c);

    fprintf('%-10s %6d %10.2f %8.2f %12.4g %10.2f %10.2f %12.4g\n', ...
        groupName{g}, numel(b), mean(b), std(b)/sqrt(numel(b)), pMeanDX(g), ...
        mean(c), std(c)/sqrt(numel(c)), pAUC(g));
end

fprintf('\nBetween-group comparison (Wilcoxon rank-sum)\n');
pairs = [1 2; 1 3; 2 3];
for k = 1:size(pairs,1)
    g1 = pairs(k,1); g2 = pairs(k,2);
    pB = ranksum(meanDX{g1}(~isnan(meanDX{g1})), meanDX{g2}(~isnan(meanDX{g2})));
    pC = ranksum(aucDX{g1}(~isnan(aucDX{g1})),   aucDX{g2}(~isnan(aucDX{g2})));
    fprintf('  %-10s vs %-10s   p (dX) = %-10.4g p (AUC) = %.4g\n', ...
        groupName{g1}, groupName{g2}, pB, pC);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  Figure
%  ------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1000 280]);

% --- Panel a ------------------------------------------------------------
ax1 = subplot(1, 3, 1); hold on;

yl = [-10 15];
patch([0 2 2 0], [yl(1) yl(1) yl(2) yl(2)], [1 0.95 0.6], ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');

for g = 1:nGroup
    fill([timePlot fliplr(timePlot)], ...
         [dxMean(g,:)+dxSEM(g,:) fliplr(dxMean(g,:)-dxSEM(g,:))], ...
         colGroup{g}, 'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
end
for g = 1:nGroup
    plot(timePlot, dxMean(g,:), '-o', 'Color', colGroup{g}, ...
         'MarkerFaceColor', colGroup{g}, 'MarkerSize', 4, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (n = %d)', groupName{g}, groupN(g)));
end

xlabel('Time from shock (s)');
ylabel('dX from baseline [+ = wall]');
text(1, yl(2)-1, 'Foot shock', 'Color', [0.85 0.2 0.2], ...
     'HorizontalAlignment', 'center', 'FontSize', 9);
xlim([timePlot(1)-1 timePlot(end)+1]); ylim(yl);
legend('Location', 'northwest', 'Box', 'on', 'FontSize', 8);
box off; set(ax1, 'TickDir', 'out');
title('a', 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

% --- Panels b and c -----------------------------------------------------
panelData  = {meanDX, aucDX};
panelP     = {pMeanDX, pAUC};
panelLabel = {'dX over 0-10 s [+ = wall]', 'AUC over 0-10 s [+ = wall]'};
panelTitle = {'b', 'c'};

for p = 1:2
    ax = subplot(1, 3, p+1); hold on;

    m = nan(1, nGroup); e = nan(1, nGroup);
    for g = 1:nGroup
        v = panelData{p}{g};
        v = v(~isnan(v));
        m(g) = mean(v);
        e(g) = std(v) / sqrt(numel(v));
    end

    for g = 1:nGroup
        bar(g, m(g), 0.6, 'FaceColor', colGroup{g}, 'EdgeColor', 'k', ...
            'LineWidth', 0.8);
    end
    errorbar(1:nGroup, m, e, 'k', 'LineStyle', 'none', 'LineWidth', 1);

    for g = 1:nGroup
        star = sigStar(panelP{p}(g));
        if ~isempty(star)
            text(g, m(g) + e(g) + 0.06*range(ylim), star, ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end

    yline(0, 'k-', 'LineWidth', 0.5);
    set(ax, 'XTick', 1:nGroup, 'XTickLabel', groupName, 'TickDir', 'out');
    ylabel(panelLabel{p});
    xlim([0.4 nGroup+0.6]);
    box off;
    title(panelTitle{p}, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
end

print(fig, 'SFig4', '-dpng', '-r300');
fprintf('Figure saved as SFig4.png\n');

%% ------------------------------------------------------------------------
%  Helper
%  ------------------------------------------------------------------------
function s = sigStar(p)
    if     isnan(p),   s = '';
    elseif p < 0.001,  s = '***';
    elseif p < 0.01,   s = '**';
    elseif p < 0.05,   s = '*';
    else,              s = 'n.s.';
    end
end
