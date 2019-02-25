%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ecModel_script
%
% - Reconstructs ecScoGEM from scoGEM
% - Make M145 and M1152 specific models
% - Simulations with proteome pool
% - Prepare proteomics data
% - Make sample-specific proteome constrained models
%
% Eduard Kerkhoven, 2019-01-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');

%% - Reconstructs ecScoGEM from scoGEM
model       = importModel([root '/ModelFiles/xml/scoGEM.xml']);
model       = setParam(model,'eq',{'BIOMASS_SCO','PROTEIN_PSEUDO'},0);
model       = setParam(model,'obj','BIOMASS_SCO_tRNA',1);

cd([ecDir '/reconstruction'])
model       = preprocessModel(model);
cd ([ecDir '/gecko/geckomat/get_enzyme_data'])
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,'streptomyces coelicolor');
cd ../change_model
ecModel    = readKcatData(model_data,kcats);
save([root '/scrap/ecScoGEM_preManualMods.mat'],'ecModel','model_data','kcats');
%load([root '/scrap/ecScoGEM_preManualMods.mat']);

cd([ecDir '/reconstruction'])
[ecModel,modifications] = manualModifications(ecModel);
save([root '/scrap/ecScoGEM_postManualMods.mat'],'ecModel','model_data','kcats');
ecModel = setParam(ecModel,'eq',{'EX_glc__D_e','EX_glu__L_e'},0); %no excretion of glutamate or glucose

%% - Make M145 and M1152 specific models
% For M1152, remove reactions involving act, red, cpk and cda clusters
% Get genes from M1152 paper (Gomez-Escribano & Bibb, 2011;
% doi:10.1111/j.1751-7915.2010.00219.x)
actGenes = cellstr(strcat('SCO', string(5074:5089)));
redGenes = cellstr(strcat('SCO', string(5879:5895)));
cpkGenes = cellstr(strcat('SCO', string(6272:6286)));
cdaGenes = cellstr(strcat('SCO', string(3212:3246)));
allGenes = [actGenes,redGenes,cpkGenes,cdaGenes];
allGenes = allGenes(contains(allGenes,model.genes));
modelM1152 = removeGenes(model,allGenes,true,true,true);
clear actGenes redGenes cpkGenes cdaGenes allGenes model

cd([ecDir '/gecko/geckomat/utilities'])
ecModel_M1152 = getSubset_ecModel(modelM1152, ecModel);

rates.M145 = dlmread([root '/ComplementaryData/growth/M145/average/M145_estimated_rates.csv'],',',1,0);
GlcUptake.M145  = -rates.M145(1,4);%
GluUptake.M145  = -rates.M145(1,5);%
gRate.M145      = rates.M145(1,3);%
ecModel=setParam(ecModel,'ub',{'EX_glc__D_e_REV',...
    'EX_glu__L_e_REV','EX_nh4_e_REV'},[GlcUptake.M145,GluUptake.M145,0]);

rates.M1152 = dlmread([root '/ComplementaryData/growth/M1152/average/M1152_estimated_rates.csv'],',',1,0);
GlcUptake.M1152  = -rates.M1152(1,4);%
GluUptake.M1152  = -rates.M1152(1,5);%
gRate.M1152      = rates.M1152(1,3);%
ecModel_M1152=setParam(ecModel_M1152,'ub',{'EX_glc__D_e_REV',...
    'EX_glu__L_e_REV','EX_nh4_e_REV'},[GlcUptake.M1152,GluUptake.M1152,0]);

load([root '/scrap/ecScoGEM_M145M1152.mat'],'ecModel*');

%% - Simulations with proteome pool
sigma   = 0.5;
Ptot    = 0.456; % As reported by Shabab et al. (1996) Microbiol.

cd([ecDir '/gecko/geckomat/limit_proteins'])

[ecModel_M145_pool, optSigma] = getConstrainedModel(ecModel,'',sigma,Ptot,gRate.M1152,modifications,'scoGEM');
sol=solveLP(ecModel_M145_pool)

[ecModel_M1152_pool, optSigma] = getConstrainedModel(ecModel_M1152,'',sigma,Ptot,gRate.M1152,modifications,'scoGEM');
sol=solveLP(ecModel_M1152_pool)

save([root '/scrap/ecPool.mat'],'ecModel*_pool')
exportModel(ecModel_M145_pool,[root,'/ModelFiles/xml/ecM145_pool.xml']);
exportModel(ecModel_M1152_pool,[root,'/ModelFiles/xml/ecM1152_pool.xml']);
end
%% - Prepare proteomics data
strain      = [string(repmat('M145',9,1));string(repmat('M1152',8,1))];
time        = string([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]);
sample      = strcat(strain,'_',time);

sigma   = 1;
Ptot    = 0.456; % As reported by Shabab et al. (1996) Microbiol.
cd([ecDir '/gecko/geckomat/limit_proteins'])
[f,~]   = measureAbundance(ecModel.enzymes);
cd(ecDir)
data    = prepareProteomicsM1152M145;
pIDs    = data.genes;
clear gRate modifications
%% - Make sample-specific proteome constrained models
cd([ecDir '/gecko/geckomat/limit_proteins'])
ecModel_M1152 = setParam(ecModel_M1152,'eq',{'EX_glc__D_e','EX_glu__L_e'},0); %no excretion of glutamate or glucose
ecModel_M1152 = setParam(ecModel_M1152,'eq',{'PSEUDO_ACCEPTOR_NAD',...
    'PSEUDO_ACCEPTOR_NADP','PSEUDO_DONOR_NADH','PSEUDO_DONOR_NADPH'},0);
ecModel = setParam(ecModel,'eq',{'EX_glc__D_e','EX_glu__L_e'},0); %no excretion of glutamate or glucose
ecModel = setParam(ecModel,'eq',{'PSEUDO_ACCEPTOR_NAD',...
    'PSEUDO_ACCEPTOR_NADP','PSEUDO_DONOR_NADH','PSEUDO_DONOR_NADPH'},0);
for i=1:length(sample)
    disp(['Generating ecModel for sample: ' sample(i)])
    protData = data.mean(:,i) + data.std(:,i);
    if i<11
        model{i} = constrainEnzymes(ecModel,Ptot,sigma,f,[],pIDs,protData);
    else
        model{i} = constrainEnzymes(ecModel_M1152,Ptot,sigma,f,[],pIDs,protData);
    end
    cd ([ecDir '\simulation'])
    [model{i}, gRate(i)] = fixFluxes(model{i},strain{i},time{i});
    if gRate(i)<0
        gRate(i)=0.0070;
    end
    cd([ecDir '/gecko/geckomat/limit_proteins'])
    [model{i},~,modifications{i}] = flexibilizeProteins(model{i},gRate(i));
    model{i} = setParam(model{i},'lb','BIOMASS_SCO_tRNA',0.99*gRate(i));
    model{i} = setParam(model{i},'ub','BIOMASS_SCO_tRNA',1.01*gRate(i));
    model{i} = setParam(model{i},'obj','ATPM',1);
end

%% All models should have the same flexibilization. First list
clear allModifications
allModifications = [modifications{1}.protein_IDs modifications{1}.modified_values];
for i=1:length(modifications)
    [Lia, Locb] = ismember(modifications{i}.protein_IDs,allModifications(:,1));
    c = {cell2mat(modifications{i}.modified_values(Lia));cell2mat(allModifications(Locb(Lia),2))};
    allModifications(Locb(Lia),2) = num2cell(max([c{:}], [], 2));
    allModifications = [allModifications;table2cell(modifications{i}(~Lia,[1,3]))];
end

for i=1:length(model)
    [Lia, Locb] = ismember(strcat(allModifications(:,1),'_exchange'),model{i}.rxnNames);
    model{i}.ub(Locb(Lia)) = cell2mat(allModifications(Lia,2));
end

for i=1:length(sample)
    model{i} = setParam(model{i},'obj','BIOMASS_SCO_tRNA',1);
    sol=solveLP(model{i},1);
    disp(['Sol = ' num2str(sol.f)])
end

save([root '/scrap/proteomeModels.mat'],'model','gRate','modifications','sample');

%% Export SBML files
for i=1:length(sample(:,1))
    exportModel(model{i},...
    [root,'/ModelFiles/xml/ec',sample{i},'.xml']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%