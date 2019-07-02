%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct_ecSco-GEM
%
% - Reconstructs ecSco-GEM from scoGEM
% - Make M145 and M1152 specific models
% - Simulations with proteome pool
%
% Eduard Kerkhoven, 2019-03-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');

%% - Reconstructs ecScoGEM from scoGEM
model      = importModel([root '/ModelFiles/xml/scoGEM.xml']);

% ec-models represent M145 and M1152 strains, these don't contain the two
% plasmids. Any SCP* gene should be removed.
scpGene    = model.genes(contains(model.genes,'SCP'));
model      = removeGenes(model,scpGene,true,true,true);

cd([ecDir '/reconstruction'])
model      = preprocessModel(model);
save([root '/scrap/Sco-GEM.mat'],'model');

cd ([ecDir '/gecko/geckomat/get_enzyme_data'])
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,'streptomyces coelicolor');
cd ../change_model
ecModel    = readKcatData(model_data,kcats);
save([root '/scrap/ecScoGEM_preManualMods.mat'],'ecModel','model_data','kcats');

cd([ecDir '/reconstruction'])
[ecModel,modifications] = manualModifications(ecModel);
save([root '/scrap/ecScoGEM_postManualMods.mat'],'ecModel','model_data','kcats');

%% - Make M145 and M1152 specific models
% For M1152, remove reactions involving act, red, cpk and cda clusters
% Get genes from M1152 paper (Gomez-Escribano & Bibb, 2011;
% doi:10.1111/j.1751-7915.2010.00219.x)
actGenes = cellstr(strcat('SCO', string(5073:5090)));
redGenes = cellstr(strcat('SCO', string(5878:5896)));
cpkGenes = cellstr(strcat('SCO', string(6270:6286)));
cdaGenes = cellstr(strcat('SCO', string(3211:3247)));
allGenes = [actGenes,redGenes,cpkGenes,cdaGenes];
allGenes = allGenes(contains(allGenes,model.genes));
modelM1152 = removeGenes(model,allGenes,true,true,true);
clear actGenes redGenes cpkGenes cdaGenes allGenes model

cd([ecDir '/gecko/geckomat/utilities'])
ecModel_M1152 = getSubset_ecModel(modelM1152, ecModel);
ecModel_M1152 = setParam(ecModel_M1152,'eq',{'EX_glc__D_e','EX_glu__L_e','EX_o2_e'},0);

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

save([root '/scrap/ecScoGEM_M145M1152.mat'],'ecModel*','gRate');

%% - Simulations with proteome pool
sigma   = 0.5;
Ptot    = 0.456; % As reported by Shabab et al. (1996) Microbiol.

cd([ecDir '/gecko/geckomat/limit_proteins'])

[ecModel_M145_pool, optSigma] = getConstrainedModel(ecModel,'',sigma,Ptot,gRate.M145,[],'scoGEM');
sol=solveLP(ecModel_M145_pool)

[ecModel_M1152_pool, optSigma] = getConstrainedModel(ecModel_M1152,'',sigma,Ptot,gRate.M1152,[],'scoGEM');
sol=solveLP(ecModel_M1152_pool)

save([root '/scrap/ecPool.mat'],'ecModel*_pool')
exportForGit(ecModel_M145_pool,'ecM145_pool',root,'xml');
exportForGit(ecModel_M1152_pool,'ecM1152_pool',root,'xml');
movefile([root '\ModelFiles\dependencies.txt'],[root '/ModelFiles/ecScoGEM_dependencies.txt'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%