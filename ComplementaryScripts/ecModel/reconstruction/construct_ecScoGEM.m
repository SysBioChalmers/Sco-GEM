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
%load([root '/scrap/ecScoGEM_postManualMods.mat'])
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

save([root '/scrap/ecScoGEM_M145M1152.mat'],'ecModel*');
load([root '/scrap/ecScoGEM_M145M1152.mat']);

%% - Simulations with proteome pool
sigma   = 0.5;
Ptot    = 0.456; % As reported by Shabab et al. (1996) Microbiol.

cd([ecDir '/gecko/geckomat/limit_proteins'])

[ecModel_M145_pool, optSigma] = getConstrainedModel(ecModel,'',sigma,Ptot,gRate.M1152,modifications,'scoGEM');
sol=solveLP(ecModel_M145_pool)

[ecModel_M1152_pool, optSigma] = getConstrainedModel(ecModel_M1152,'',sigma,Ptot,gRate.M1152,modifications,'scoGEM');
sol=solveLP(ecModel_M1152_pool)

% conditions  = num2cell(sample(:,[1,3]),1);
% conditions  = strcat(conditions{:,1},{'_'},conditions{:,2});
% T = topUsedEnzymes(transpose(sol_pool),model_pool(1),conditions,'scoGEM',false,10);
% writetable(T,['../../../../../ComplementaryData/ecmodel/ecScoGEM_pool_topUsedEnzymes.txt'])

save([root '/scrap/ecPool.mat'],'ecModel*_pool')
load([root '/scrap/ecPool.mat'])

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
ecModel = setParam(ecModel,'eq',{'EX_glc__D_e','EX_glu__L_e'},0); %no excretion of glutamate or glucose
for i=1:length(sample)
    disp(['Generating ecModel for sample: ' sample(i)])
    protData = data.mean(:,i) + data.std(:,i);
    model.(sample(i)) = constrainEnzymes(ecModel,Ptot,sigma,f,[],pIDs,protData);
    cd ([ecDir '\simulation'])
    [model.(sample(i)), gRate(i)] = fixFluxes(model.(sample(i)),strain{i},time{i});
    if gRate(i)<0
       gRate(i)=0.0070;
    end    
   cd([ecDir '/gecko/geckomat/limit_proteins'])
   [model.(sample(i)),~,modifications.(sample(i))] = flexibilizeProteins(model.(sample(i)),gRate(i));
end

save([root '/scrap/proteomeModels.mat'],'model','gRate','modifications');

%% Estimate ATP maintenance
modelCell=struct2cell(model);
for i=1:length(gRate)
    modelCell{i} = setParam(modelCell{i},'obj','BIOMASS_SCO_tRNA',1);
    sol(i) = solveLP(modelCell{i});
    tmp = setParam(modelCell{i},'eq','BIOMASS_SCO_tRNA',-sol(i).f);
    tmp = setParam(tmp,'obj','ATPM',1);
    tmp2 = solveLP(tmp);
    tmp = setParam(tmp,'obj','BIOMASS_SCO_tRNA',1);
    modelCell{i}= setParam(tmp,'eq','ATPM',-0.999*tmp2.f);
end

%% Export SBML files
for i=1:length(sample(:,1))
    exportModel(modelCell{i},...
    [root,'/ModelFiles/xml/ec',sample{i},'.xml']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%