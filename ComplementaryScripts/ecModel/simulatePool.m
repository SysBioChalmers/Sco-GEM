%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulatePool(ecModel)
%
% Eduard Kerkhoven, 2018-12-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate model
load('..\..\ModelFiles\mat\ECscoGEM.mat')
data        = prepareProteomicsM1152M145;
sol         = solveLP(ecModel)

cd gecko/geckomat/limit_proteins
Ptot        = 0.412;
sigma       = 0.42;
[f,~]       = measureAbundance(ecModel.enzymes);
ecModel_pool = constrainPool(ecModel,true(length(ecModel.enzymes),1),f*sigma*Ptot);
ecModel_pool = setParam(ecModel_pool,'lb','ATPM',0);

cd ../../../
exportForGit(ecModel_pool,'ecScoGEM_pool','../../',{'xml','txt','yml','xlsx'})

%% Simulate with experimental data
sample      = [cellstr(repmat('M145',9,1));cellstr(repmat('M1152',8,1))];
sample(:,2) = num2cell([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]);
sample(:,3) = cellstr(num2str([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]));

for i=1:length(sample(:,1))%[1:3,10:12];%
    model_pool(i)   = simulateCondition(ecModel_pool,sample(i,1),sample{i,2});
    sol = solveLP(model_pool(i),1);
    sol_pool(i,:) = sol.x;
end

%% Plot which are the top used enzymes
cd gecko/geckomat/kcat_sensitivity_analysis

conditions  = num2cell(sample(:,[1,3]),1);
conditions  = strcat(conditions{:,1},{'_'},conditions{:,2});
topUsedEnzymes(transpose(sol_pool),model_pool(1),conditions,'scoGEM',true,1000);