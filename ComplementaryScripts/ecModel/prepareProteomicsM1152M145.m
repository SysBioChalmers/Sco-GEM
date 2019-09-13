%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = prepareProteomicsM1152M145
%
% Eduard Kerkhoven, 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = prepareProteomicsM1152M145
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

load('../ecModel/gecko/databases/ProtDatabase.mat')
[~, ib]=ismember(data.genes,kegg(:,3));
data.MW     = cell2mat(kegg(ib,5))./1000; % kDa

%% Normalize for size, and convert to mmol/g protein
data.norm   = data.top3 .* data.MW; % correct for each protein's kDa
data.norm   = data.norm ./ sum(data.norm,1,'omitnan'); % weight ratio (g/g)
data.norm   = data.norm ./ data.MW; % mmol / g protein

%data.norm   = data.norm .* 0.0002; % injected 0.2 ug, new unit: mg/0.2 ug
%data.norm   = data.norm ./ (data.MW*1000); % convert to mmol/0.2 ug

%% Assume 0.4 g protein / gDCW
data.norm   = data.norm .* 0.456; % As reported by Shabab et al. (1996)
                % Microbiol. doi:10.1099/13500872-142-8-1927 relatively
                % stabile among growth rates, highest reported content is
                % used to prevent overconstraining the model.

%% Average over 3 replicates, M145 has 9 timepoints, M1152 has 8 timepoints
for i=1:9
    data.mean(:,i)  = mean(data.norm(:,[i,i+9,i+18]),2,'omitnan');
    data.std(:,i)   = std(data.norm(:,[i,i+9,i+18]),0,2,'omitnan');
end
for i=10:17
    j=i+18;
    data.mean(:,i)  = mean(data.norm(:,[j,j+8,j+16]),2,'omitnan'); 
    data.std(:,i)   = std(data.norm(:,[j,j+8,j+16]),0,2,'omitnan');
end

data.meanSample     = unique(regexprep(data.sample,'_F...',''),'stable');
end