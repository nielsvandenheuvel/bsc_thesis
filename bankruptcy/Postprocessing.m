%% Helper Functions
iff = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
curly = @(x, varargin) x{varargin{:}};

%% Load Survival Data
hbartSurv = readtable("./Code/bankruptcy/heter_surv.csv");
bartSurv = readtable("./Code/bankruptcy/homo_surv.csv");
dataArray = readtable("Bankruptcy-Bart.csv");

dataArray = dataArray(~sum(isnan(dataArray{:, :}), 2) > 0, :);
dataArray = dataArray(dataArray.t < 2021, :);

%% Format Data
gvkey = unique(hbartSurv.gvkey);
gvkeyBankrupt = unique(dataArray.gvkey(dataArray.y == 1));
gvkeyLive = setdiff(gvkey, gvkeyBankrupt);
numFirms = length(gvkey);

data(numFirms) = struct();
for i = 1:numFirms
    dataRows = find(hbartSurv{:, "gvkey"} == gvkey(i));
    
    data(i).gvkey = gvkey(i);
    data(i).status = ismember(gvkey(i), gvkeyBankrupt);
    data(i).heter = hbartSurv(dataRows, 3:6);
    data(i).homo = bartSurv(dataRows, 3:6);
    
    data(i).heter.t = data(i).heter.t + 2000;
    data(i).homo.t = data(i).homo.t + 2000;
    
    data(i).heter.y = dataArray.y(find(dataArray.gvkey == gvkey(i)));
    data(i).homo.y = dataArray.y(find(dataArray.gvkey == gvkey(i)));
end

%% Test for differences in means for covariates
test = struct();
ratios = ["pe_inc", "bm", "cash_lt", "ps", "debt_at", "lt_debt", "roe", "short_debt"];
for i = 1:8
    test(i).ratio = ratios(i);
    test(i).mean = [mean(unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}), 'omitnan'), mean(unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}), 'omitnan')];
    test(i).sd = [std(unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}), 'omitnan'), std(unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}), 'omitnan')];
    [test(i).h, test(i).p, test(i).ci, test(i).stats] = ttest2(unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}), unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}), 'VarType', 'unequal');
end

%% Plot boxplots
tiledlayout(2,4)
ylimits = [-30 50; 0 3; 0 2; 0 6; 0 1; 0 1; -1 1; 0 1];
titles = ["Price/Equity", "Book/Market", "Cash Balance/Total Liabilities", "Price/Sales", "Total Debt/Total Assets", "Long-term Debt/Total Liabilities", "Return on Equity", "Short-term Debt/Total Debt"];
for i = 1:length(ratios)
    nexttile
    g = [repmat(1,length(unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio})),1); repmat(2,length(unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio})),1)];
    x = [unique(dataArray{ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio}); unique(dataArray{~ismember(dataArray.gvkey, gvkeyBankrupt), test(i).ratio})];
    boxplot(x, g, 'Symbol', '')
    ylim(ylimits(i, :))
    title(strcat(titles(i), num2str(test(i).p)))
end

%% Put data in a time-series structure
panel(length(data(1).heter.t)) = struct();
for t = 1:length(panel)
    
    panel(t).t = data(1).heter.t(t);
    
    if t == 1
        status = arrayfun(@(x) iff(ismember(panel(t).t, x.heter.t), x.heter.y(x.heter.t == panel(t).t), true, @() -1), data)';
    else
        status = arrayfun(@(x) iff(ismember(panel(t).t, x.heter.t), x.heter.y(x.heter.t == panel(t).t), true, @() -1), data)';
        status = arrayfun(@(i) iff(panel(t - 1).data.status(i), -1, true, status(i)), 1:length(status))';
    end
    gvkey = [data.gvkey]';
    
    heter_p = arrayfun(@(x) iff(ismember(panel(t).t, x.heter.t), x.heter.p(x.heter.t == panel(t).t), true, @() NaN), data)';
    heter_S = arrayfun(@(x) iff(ismember(panel(t).t, x.heter.t), x.heter.S(x.heter.t == panel(t).t), true, @() NaN), data)';
    
    homo_p = arrayfun(@(x) iff(ismember(panel(t).t, x.homo.t), x.homo.p(x.homo.t == panel(t).t), true, @() NaN), data)';
    homo_S = arrayfun(@(x) iff(ismember(panel(t).t, x.homo.t), x.homo.S(x.homo.t == panel(t).t), true, @() NaN), data)';
    
    panel(t).data = array2table([gvkey, status, heter_p, heter_S, homo_p, homo_S],...
                                "VariableNames", ["gvkey", "status", "heter_p", "heter_S", "homo_p", "homo_S"]);
end

%% Kaplan-Meier Estimator
kaplanMeier = table();
kaplanMeier.t = data(1).heter.t;
kaplanMeier.p = arrayfun(@(t) sum(t.data.status == 1)/sum(t.data.status >= 0), panel)';
for t = 1:length(kaplanMeier.t)
   kaplanMeier.S(t) = prod(1 - kaplanMeier.p(1:t));
end

%% ROC Curve
numBankrupt = sum(arrayfun(@(s) sum(panel(s).data.status == 1), 1:length(panel)));
heter_fracFailed = zeros(101, length(panel));
for t = 1:length(panel)
    % Get the probabilities for the firms that are still in the sample.
    [~, I] = sort(panel(t).data.heter_p);
    binData = panel(t).data(I, :);
    binData = binData(binData.status >= 0, :);
    % Create bins
    [idBin, bins] = discretize(binData.heter_p, 100);
    for i = 1:length(bins)
        heter_fracFailed(i, t) = sum(binData.status(idBin <= i))/numBankrupt;
    end
end

homo_fracFailed = zeros(101, length(panel));
for t = 1:length(panel)
    % Get the probabilities for the firms that are still in the sample.
    [~, I] = sort(panel(t).data.homo_p);
    binData = panel(t).data(I, :);
    binData = binData(binData.status >= 0, :);
    % Create bins
    [idBin, bins] = discretize(binData.homo_p, 100);
    for i = 1:length(bins)
        homo_fracFailed(i, t) = sum(binData.status(idBin <= i))/numBankrupt;
    end
end

figure
plot(0:0.01:1, sum(heter_fracFailed, 2), 'k')
hold on
plot(0:0.01:1, sum(homo_fracFailed, 2), '--k')
plot(0:0.01:1,0:0.01:1, ':k')
hold off
ax = gca;
ax.FontSize = 12;
xlabel("% of Firms", 'FontSize', 12)
ylabel("% of Failed Firms", 'FontSize', 12)
grid on
% Are hazard models superior to traditional bankruptcy prediction approaches? A comprehensive test

%% Survival Function
% Bankrupt vs non-bankrupt
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
plot(arrayfun(@(t) median(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 'omitnan'), panel), 'o-k', 'MarkerFaceColor', 'w')
hold on
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.25), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.75), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.05), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_S(ismember(t.data.gvkey, gvkeyBankrupt)), 0.95), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold off
grid on
ylim([0 1])
xlim([1 21])
xticks(1:2:21)
xticklabels(2000:2:2020)
ax = gca;
ax.FontSize = 12;
title("Bankrupt", 'FontSize', 12)
xlabel("Year", 'FontSize', 12)
ylabel("Survival Probability", 'FontSize', 12)

nexttile
plot(arrayfun(@(t) median(t.data.heter_S(~ismember(t.data.gvkey, gvkeyBankrupt)), 'omitnan'), panel), 'o-k', 'MarkerFaceColor', 'w')
hold on
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_S(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.25), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_S(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.75), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_S(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.05), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_S(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.95), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold off
grid on
ylim([0 1])
xlim([1 21])
xticks(1:2:21)
xticklabels(2000:2:2020)
ax = gca;
ax.FontSize = 12;
title("Viable", 'FontSize', 12)
xlabel("Year", 'FontSize', 12)

%% Brier Scores
brierScore.heter = arrayfun(@(t) mean((t.data.heter_p(t.data.status >= 0) - t.data.status(t.data.status >= 0)).^2), panel);
brierScore.homo = arrayfun(@(t) mean((t.data.homo_p(t.data.status >= 0) - t.data.status(t.data.status >= 0)).^2), panel);

figure
b = bar([brierScore.heter; brierScore.homo]');
hold on 
plot(arrayfun(@(t) sum(t.data.status(t.data.status >= 0))/sum(t.data.status >= 0), panel), '-.k')
hold off
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0.5 0.5];
xticks(1:2:21)
xticklabels(2000:2:2020)
yticks(0:0.025:0.2)
grid on
xlabel("Year")
ylabel("Brier score")


%% Miscellaneous Figures
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
plot(arrayfun(@(t) median(t.data.heter_p(ismember(t.data.gvkey, gvkeyBankrupt)), 'omitnan'), panel), 'o-k', 'MarkerFaceColor', 'w')
hold on
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_p(ismember(t.data.gvkey, gvkeyBankrupt)), 0.25), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_p(ismember(t.data.gvkey, gvkeyBankrupt)), 0.75), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_p(ismember(t.data.gvkey, gvkeyBankrupt)), 0.05), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_p(ismember(t.data.gvkey, gvkeyBankrupt)), 0.95), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold off
grid on
ylim([0 0.3])
xlim([1 21])
xticks(1:2:21)
xticklabels(2000:2:2020)
ax = gca;
ax.FontSize = 12;
title("Bankrupt", 'FontSize', 12)
xlabel("Year", 'FontSize', 12)
ylabel("Survival Probability", 'FontSize', 12)

nexttile
plot(arrayfun(@(t) median(t.data.heter_p(~ismember(t.data.gvkey, gvkeyBankrupt)), 'omitnan'), panel), 'o-k', 'MarkerFaceColor', 'w')
hold on
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_p(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.25), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_p(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.75), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
fill([1:21, fliplr(1:21)], [arrayfun(@(t) quantile(t.data.heter_p(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.05), panel), fliplr(arrayfun(@(t) quantile(t.data.heter_p(~ismember(t.data.gvkey, gvkeyBankrupt)), 0.95), panel))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold off
grid on
ylim([0 0.3])
xlim([1 21])
xticks(1:2:21)
xticklabels(2000:2:2020)
ax = gca;
ax.FontSize = 12;
title("Viable", 'FontSize', 12)
xlabel("Year", 'FontSize', 12)