%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = reorderBiomass(model)
%
% Split biomass equation into smaller sub-equations describing certain
% classes of biomass components, allowing them to be scaled individually.
% 
% 
% Eduard Kerkhoven, 2018-10-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = reorderBiomass(model)

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
BM.S(contains(BM.met,'misc_c')) = 1;

%% Write new biomass coefficients in file
writetable(struct2table(BM),'../../ComplementaryData/biomass/biomass_scaled.txt','Delimiter','\t');

psIdx=find(contains(BM.name,'pseudometabolite'));
[~,psUniq,~]=unique(BM.name(psIdx));
psIdx=psIdx(psUniq);

metsToAdd.metNames=BM.name(psIdx);
metsToAdd.mets=BM.met(psIdx);
metsToAdd.compartments='c';
model=addMets(model,metsToAdd);

rxnsToAdd.rxnNames={'protein pseudoreaction','lipid pseudoreaction',...
    'DNA pseudoreaction','RNA pseudoreaction','carbohydrate pseudoreaction',...
    'misc pseudoreaction','biomass pseudoreaction','cell wall pseudoreaction'};
rxnsToAdd.rxns={'BM_prot','BM_lipid','BM_DNA','BM_RNA','BM_carb',...
    'BM_misc','BM_growth','BM_cellwall'};
     
rxnsToAdd.mets={...
    BM.met(strcmp(BM.rxn,'P')), BM.met(strcmp(BM.rxn,'L')),...
    BM.met(strcmp(BM.rxn,'D')), BM.met(strcmp(BM.rxn,'R')),...
    BM.met(strcmp(BM.rxn,'C')), BM.met(strcmp(BM.rxn,'M')),...
    BM.met(strcmp(BM.rxn,'B')), BM.met(strcmp(BM.rxn,'W'))};
rxnsToAdd.stoichCoeffs={...
    BM.S(strcmp(BM.rxn,'P')), BM.S(strcmp(BM.rxn,'L')),...
    BM.S(strcmp(BM.rxn,'D')), BM.S(strcmp(BM.rxn,'R')),...
    BM.S(strcmp(BM.rxn,'C')), BM.S(strcmp(BM.rxn,'M')),...
    BM.S(strcmp(BM.rxn,'B')), BM.S(strcmp(BM.rxn,'W'))};   
rxnsToAdd.lb=repmat(0,8,1);
model=addRxns(model,rxnsToAdd);
model=setParam(model,'eq','BIOMASS_SCO',0);
end