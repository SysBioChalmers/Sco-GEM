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
Ptot        = 0.456; % As reported by Shabab et al. (1996) Microbiol.
                     % doi:10.1099/13500872-142-8-1927 relatively stabile
                     % among growth rates, highest reported content is used
                     % to prevent overconstraining the model.
sigma       = 0.5;
[f,~]       = measureAbundance(ecModel.enzymes);
rates = dlmread('../../../../../ComplementaryData/growth/M145_estimated_rates.csv',';',1,0);
gRate      = rates(1,3);%

[ecModel_pool, optSigma] = getConstrainedModel(ecModel,'',sigma,Ptot,gRate,modifications,'scoGEM');

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
T = topUsedEnzymes(transpose(sol_pool),model_pool(1),conditions,'scoGEM',false,10);
writetable(T,['../../../../../ComplementaryData/ecmodel/' name '_topUsedEnzymes.txt'])