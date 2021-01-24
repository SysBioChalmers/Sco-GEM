%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and organize prosthetic group data
dat = readtable('../../ComplementaryData/biomass/prosthetic_groups_uniProt.txt',...
    'ReadVariableNames',0);

pGroups.genes   = table2array(dat(2:end,1));
pGroups.mets    = transpose(table2array(dat(1,2:end)));
pGroups.S       = str2double(table2array(dat(2:end,2:end)));
clear dat

%% Distribute "metal" and discard dipyrromethane
allMetals           = contains(pGroups.mets,{'cobalt2_c','cu2_c','fe2_c',...
    'zn2_c','ni2_c','ca2_c','k_c','mg2_c','mn2_c'});
allMetals           = allMetals./sum(allMetals);

metalIdx            = strcmp(pGroups.mets,'metal');
metal               = find(pGroups.S(:,metalIdx));
pGroups.S(metal,:)  = pGroups.S(metal,:) + ...
                        transpose(repmat(allMetals,1,length(metal)));
pGroups.S(:,metalIdx)   = [];
pGroups.mets(metalIdx)  = [];

% Dipyrromethane is generated by the enzyme itself from its substrate
dpmIdx              = strcmp(pGroups.mets,'dpm');
pGroups.S(:,dpmIdx)	= [];
pGroups.mets(dpmIdx)= [];

%% Load and organize proteomics data
cd ../ecModel/
[prot.data, ~, prot.genes] = prepareProteomicsM1152M145;
prot.mean           = mean(prot.data,2,'omitnan');
[Lia, Locb]         = ismember(pGroups.genes,prot.genes);
pGroups.mmol        = zeros(length(pGroups.genes),1)
pGroups.mmol(Lia)   = prot.mean(Locb(Lia));

clear dat
cd ../consensusModel/

%% Calculate abundances
pGroups.coeffs = transpose(pGroups.S) * pGroups.mmol;
out = pGroups.mets;
out(:,2) = cellstr(num2str(pGroups.coeffs,'%s'));
out = table(out(:,1),out(:,2),'VariableNames',{'mets';'coeffs'});
writetable(out,'../../ComplementaryData/biomass/prosthetic_groups_mets.txt','Delimiter','\t');

%% Refit biomass components
 
fID          = fopen('../../ComplementaryData/biomass/standard_biomass.txt'); 
data_cell    = textscan(fID,'%s %s %f %f %s','delimiter','\t','HeaderLines',1); 
fclose(fID); 
 
BM.met=data_cell{1}; 
BM.rxn=data_cell{2}; 
BM.MW=data_cell{3}; 
BM.S=data_cell{4}; 
BM.name=data_cell{5}; 
 
fID          = fopen('../../ComplementaryData/biomass/prosthetic_groups_mets.txt'); 
data_cell    = textscan(fID,'%s %f','delimiter','\t','HeaderLines',1); 
fclose(fID); 
 
PG.mets=data_cell{1}; 
PG.S=data_cell{2}; 
 
% Remove metals and NAD 
nonMetals   = ~contains(PG.mets,{'cobalt2_c','cu2_c','fe2_c',... 
    'zn2_c','ni2_c','ca2_c','k_c','mg2_c','mn2_c','nad_c'}); 
 
PG.mets(~nonMetals) = []; 
PG.S(~nonMetals)    = []; 
PG.S(isnan(PG.S)) 	= 0; 
% Replace coefficients of measure co-factors 
miscIdx = find(strcmp(BM.rxn,'M')); 
[A,nonMetalIdx]     = ismember(BM.met(miscIdx),PG.mets); 
BM.S(miscIdx(A)) = -PG.S(nonMetalIdx(A)); 
 
% Adjust remaining misc metabolites to ensure that the whole biomass adds 
% up to 1 g/gDCW 
weightDiff = 1000 + sum(BM.S.*BM.MW,'omitnan'); 
 
% This is the difference that needs to be compensated for, by "misc" 
% metabolites that were not defined as co-factors 
sumMisc = -sum(BM.S(miscIdx(~A)).*BM.MW(miscIdx(~A)),'omitnan'); 
ratio   = (weightDiff+sumMisc)/sumMisc; 
BM.S(miscIdx(~A)) = BM.S(miscIdx(~A)).*ratio; 
BM.S(strcmp(BM.met,'misc_c') & strcmp(BM.rxn,'M')) = 1; 
 
%% Write new biomass coefficients in file 
writetable(struct2table(BM),'../../ComplementaryData/biomass/biomass_scaled.txt','Delimiter','\t'); 