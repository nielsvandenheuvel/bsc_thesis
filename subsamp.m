% Generates a subsample to reproduce the preprocessing code.

data = readtable("bankruptcy/WRDS_Financial_Ratios.csv");
gvkey = unique(data.gvkey);
gvkeySub = randsample(gvkey, 1200);
data(~ismember(data.gvkey, gvkeySub), :) = [];

% writetable(data, "subsamp.csv");