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

fID          = fopen('../../ComplementaryData/biomass/standard_biomass.tab');
data_cell    = textscan(fID,'%s %s %f %f %s','delimiter','\t','HeaderLines',1);
fclose(fID);

BM.met=data_cell{1};
BM.rxn=data_cell{2};
BM.MW=data_cell{3};
BM.S=data_cell{4};
BM.name=data_cell{5};

% Adjust "misc" components to ensure that the whole biomass adds up to
% 1 g/gDCW
BM.S(strcmp(BM.rxn,'M'))=BM.S(strcmp(BM.rxn,'M'))*0.75536; 
BM.S(strcmp(BM.met,'misc_c') & strcmp(BM.rxn,'M'))=1; 


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