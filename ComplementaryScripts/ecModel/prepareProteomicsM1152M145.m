%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = prepareProteomicsM1152M145
%
% Eduard Kerkhoven, 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
dat         = readtable('../../ComplementaryData/data/proteomics_M1152M145_2018.csv');
% Pseudogene SCO7079: only in 5 samples, while no MW provided by KEGG. Just
% remove this data
idx         = find(strcmp(dat.x___entry,'SCO7079'));
dat(idx,:)  = [];
% Reorganize remaining data
data.top3   = table2array(dat(:,2:52));
data.IDs    = transpose(dat.Properties.VariableNames);
data.IDs(1) = [];
data.genes  = table2array(dat(:,1));
data.genes  = regexprep(data.genes,'(SCO\d{4}).+','$1'); % Remove trailing characters, not sure why these were included for these genes



dat         = readtable('../../ComplementaryData/data/proteomics_M1152M145_2018_ID.csv');
data.sample = strcat(dat.Strain,'_',dat.FermenterID,'_',dat.SamplingID);

load('gecko/databases/ProtDatabaseAllUniprot.mat')
[~, ib]=ismember(data.genes,kegg(:,3));
data.MW     = cell2mat(kegg(ib,5))./1000; % kDa

%% Normalize for size
data.norm   = data.top3 .* data.MW; % correct for each protein's kDa
data.norm   = data.norm ./ sum(data.norm,1,'omitnan'); % weight ratio (g/g)
data.norm   = data.norm .* 0.2; % injected protein: 0.2 ug
data.norm   = data.top3./data.MW;
%% Total count
data.size   = sum(data.norm,1,'omitnan');