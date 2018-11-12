%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ecModel = simulateProteomics(ecModel,data,sample)
%
% Eduard Kerkhoven, 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('gecko\models\ecScoGEM_pool.mat')

[data, sample, genes] = prepareProteomicsM1152M145;

cd gecko/geckomat/limit_proteins
[f,~] = measureAbundance(ecModel.enzymes,genes,data(:,1));


ecModel_pool = constrainPool(ecModel,true(length(ecModel.enzymes),1),sigma*f);



[ecModel_prot,~,~] = constrainEnzymes(ecModel,Ptot,sigma,f);

% model = constrainEnzymes(model,Ptot,sigma,f,GAM,pIDs,data,gRate,GlucUptake)


gRate       = 0.03427;
GlucUptake  = 0.51
ecModel_prot = constrainEnzymes(ecModel,0.4,0.4,f,70,genes,data(:,1),gRate,GlucUptake);




ecModel_batch = modifyKcats(ecModel_batch,obj_Val,modifications,name);
