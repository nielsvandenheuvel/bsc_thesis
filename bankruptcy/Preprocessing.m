%% Preamble
% This script preprocesses the raw data from the bankruptcy dataset and the
% finanical ratios. It starts of by formatting the data in struct and by
% iterating some criteria refines the data. Note that this file is a toy
% example, as the data used in here is only a small subsample of the whole
% data used in the analysis. The initial script wrote the entire dataset
% back to bankrutpcy/Bankruptcy-Bart.csv, which is used in the analysis.

%% Helper functions
iff = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
paren = @(x, varargin) x(varargin{:});

%% Load bankruptcy data
bankruptcyData = readtable("Florida-UCLA-LoPucki Bankruptcy Research Database 1-12-2023.csv");

% Drop variables
bankruptcyData(:, 1:71) = [];
bankruptcyData(:, 73:128) = [];
bankruptcyData(:, 130:end) = [];

% Coarsen DateFiled to YearFiled
bankruptcyData.YearFiled = year(bankruptcyData.DateFiled);
bankruptcyData.DateFiled = [];

%% Load (subsample) financial ratios data
ratiosData = readtable("subsamp.csv");

%% Format data
% Select variables
allVariables = ratiosData.Properties.VariableNames;
dropVariables = ["permno", "adate", "qdate", "TICKER", "cusip", "divyield"];
selectedCovariates = setdiff(setdiff(ratiosData.Properties.VariableNames, dropVariables), ["gvkey", "public_date"]);
ratiosData(:, dropVariables) = [];

% Delete observations before 2000
ratiosData(ratiosData{:, "public_date"} < "2000-01-01", :) = [];

% Coarsen data to years
ratiosData.public_date = year(ratiosData.public_date);

% Obtain keys
gvkey = unique(ratiosData{:, "gvkey"});
gvkeyBankrupt = intersect(bankruptcyData.GvkeyBefore(bankruptcyData.YearFiled >= 2000), gvkey);
gvkeyLive = setdiff(gvkey, gvkeyBankrupt);
rng(0);
gvkeySample = [randsample(gvkeyLive, 500); gvkeyBankrupt];

% Time plots for selected predictors.
ratios = ["pe_inc", "bm", "cash_lt", "ps", "debt_at", "lt_debt", "roe", "short_debt"];
titles = ["Price/Equity", "Book/Market", "Cash Balance/Total Liabilities", "Price/Sales", "Total Debt/Total Assets", "Long-term Debt/Total Liabilities", "Return on Equity", "Short-term Debt/Total Debt"];
tiledlayout(4,2, 'TileSpacing', 'compact')
for i = 1:length(ratios)
    nexttile
    hold on
    plot(arrayfun(@(j) mean(ratiosData{ismember(ratiosData.gvkey, gvkeyBankrupt) & ratiosData.public_date == j, ratios(i)}, 'omitnan'), 2000:2020), 'k')
    plot(arrayfun(@(j) mean(ratiosData{~ismember(ratiosData.gvkey, gvkeyBankrupt) & ratiosData.public_date == j, ratios(i)}, 'omitnan'), 2000:2020), '--k')
    hold off
    box on
    grid on
    xlabel('Year')
    xlim([1 21])
    xticks(1:4:21)
    xticklabels(2000:4:2020)
    ax = gca;
    ax.FontSize = 12;
    title(titles(i), 'Fontsize', 12)
end


numFirms = length(gvkeySample);
data(numFirms) = struct();
for i = 1:numFirms
    
    % Set rows for firm
    ratiosDataRows = find(ratiosData{:, "gvkey"} == gvkeySample(i));
    
    % Set GVKEY and status
    data(i).gvkey = gvkeySample(i);
    data(i).status = min(sum(bankruptcyData.GvkeyBefore == gvkeySample(i)), 1);
    
    % Set time vector
    T = unique(ratiosData.public_date(ratiosDataRows));
    
    % Add bankruptcy status
    y = zeros(length(T), 1);
    if data(i).status == 1
        y(length(T)) = 1;
    end
    
    % Add ratios
    ratios = NaN(length(T), length(selectedCovariates));
    for t = 1:length(T)
        if sum(ratiosData.public_date(ratiosDataRows) == T(t)) >= 1
            ratios(t, :) = mean(ratiosData{ratiosDataRows, selectedCovariates}, 1, 'omitnan');
        end
    end
    data(i).table = array2table([T, y, ratios], "VariableNames", ["t", "y", selectedCovariates]);
    
    if mod(i, 100) == 0
        disp(i);
    end
end

%% Refine data
% Only get rows that have more than 9 years of data or 4 years, if
% bankrupt.
filter = arrayfun(@(x) size(x.table, 1) >= 10 | (x.status == 1 & size(x.table, 1) >= 4), data);
filteredData = data(filter);

%% Missing data
% Compute the missing fraction of observations for each variable.
fracMissing = array2table(arrayfun(@(x) sum(arrayfun(@(i) sum(isnan(i.table{:, x})), filteredData)), 1:70)/sum(arrayfun(@(i) size(i.table, 1), filteredData)));
fracMissing.Properties.VariableNames = filteredData(1).table.Properties.VariableNames;

% These variables are often 0 and I expect them to carry little explanatory
% power.
dropMissingVariables = ["rd_sale", "staff_sale", "adv_sale"];

% Add variables that have more than 4% missing data
dropMissingVariables = [dropMissingVariables, fracMissing.Properties.VariableNames(fracMissing{:, :} > 0.04)];

% Get the bankrupt firms that have missing data at the year of bankruptcy
missingBankruptFirms = arrayfun(@(i) iff(i.status == 1, sum(isnan(i.table{size(i.table, 1), ~ismember(i.table.Properties.VariableNames, dropMissingVariables)})), true, 0), filteredData) > 0;
filteredData(missingBankruptFirms) = [];

% Recalculate missing values
variableLocations = find(~ismember(fracMissing.Properties.VariableNames, dropMissingVariables));
fracMissing = array2table(arrayfun(@(x) sum(arrayfun(@(i) sum(isnan(i.table{:, x})), filteredData)), variableLocations)/sum(arrayfun(@(i) size(i.table, 1), filteredData)));
fracMissing.Properties.VariableNames = filteredData(1).table.Properties.VariableNames(variableLocations);

% The rest can be dealt with by complete.cases in R.

%% Split into estimation and prediction sample
rng(0);
iTrain = randsample(1:length(filteredData), 0.5*length(filteredData));
iPred = setdiff(1:length(filteredData), iTrain);

for i = 1:length(filteredData)
    filteredData(i).train = ismember(i, iTrain);
    % If firm is in prediction sample, add prediction data
    if ~filteredData(i).train
        T = (filteredData(i).table.t(1):2022)';
        covData = repmat(filteredData(i).table{1, 3:end}, length(T), 1);
        filteredData(i).pred = array2table([T, covData], "VariableNames", ["t", selectedCovariates]);
    end
end

%% Convert structure back to array
dataArray = array2table(NaN(sum(arrayfun(@(x) size(x.table, 1), filteredData)), 73));
dataArray.Properties.VariableNames = ["gvkey", "train", "status", filteredData(1).table.Properties.VariableNames];
currentRow = 1;
endRow = 0;
for i = 1:length(filteredData)
    if filteredData(i).train
        endRow = endRow + size(filteredData(i).table, 1);
        dataArray(currentRow:endRow, filteredData(1).table.Properties.VariableNames) = filteredData(i).table;
    else
        endRow = endRow + size(filteredData(i).pred, 1);
        dataArray(currentRow:endRow, "y") = array2table([filteredData(i).table.y; ones((endRow - currentRow + 1) - length(filteredData(i).table.y), 1)]);
        dataArray(currentRow:endRow, filteredData(1).pred.Properties.VariableNames) = filteredData(i).pred;
    end
    dataArray.gvkey(currentRow:endRow) = filteredData(i).gvkey;
    dataArray.train(currentRow:endRow) = repmat(filteredData(i).train, endRow - currentRow + 1, 1);
    dataArray.status(currentRow:endRow) = filteredData(i).status;
    currentRow = endRow + 1;

    if mod(i, 1000) == 0
        disp(i)
    end
end

dataArray = dataArray(:, ~ismember(dataArray.Properties.VariableNames, dropMissingVariables));