%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ecModel = simulateProteomics(ecModel,data,sample)
%
% Eduard Kerkhoven, 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('..\..\ModelFiles\mat\ECscoGEM.mat')
data        = prepareProteomicsM1152M145;
sol         = solveLP(ecModel);
printFluxes(ecModel,sol.x)

rates       = dlmread('../../ComplementaryData/growth/M145_estimated_rates.csv',';',1,0);
gRate       = rates(1,3);%
cd gecko/geckomat/limit_proteins
Ptot        = 0.412;
sigma       = 0.42;
[f,~]       = measureAbundance(ecModel.enzymes);
GAM         = [];
pIDs        = data.genes;
ecModel_pool = constrainPool(ecModel,true(length(ecModel.enzymes),1),f*sigma*Ptot);

sample      = [cellstr(repmat('M145',9,1));cellstr(repmat('M1152',8,1))];
sample(:,2) = num2cell([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]);
sample(:,3) = cellstr(num2str([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]));

ecModel_pool=setParam(ecModel_pool,'lb','ATPM',0);
cd ../../../
for i=1:length(sample(:,1))%[1:3,10:12];%
    model_pool(i)   = simulateCondition(ecModel_pool,sample(i,1),sample{i,2});
    sol = solveLP(model_pool(i),1);
    sol_pool(i,:) = sol.x;
end

%% Plot which are the top used enzymes
cd gecko/geckomat/kcat_sensitivity_analysis
conditions={'M145_21','M145_29','M145_33','M145_37','M145_41','M145_45',...
    'M145_49','M145_53','M145_57','M1152_33','M1152_41','M1152_45',...
    'M1152_49','M1152_53','M1152_57','M1152_61','M1152_65'};
topUsedEnzymes(transpose(sol_pool),model_pool(1),conditions,'scoGEM');
