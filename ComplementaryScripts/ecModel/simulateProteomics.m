%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ecModel = simulateProteomics(ecModel,data,sample)
%
% Eduard Kerkhoven, 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('..\..\ModelFiles\mat\ecScoGEM.mat')
data        = prepareProteomicsM1152M145;
sol         = solveLP(ecModel)

rates       = dlmread('../../ComplementaryData/growth/M145_estimated_rates.csv',';',1,0);
gRate       = rates(1,3);%
cd gecko/geckomat/limit_proteins
Ptot        = 0.456; % As reported by Shabab et al. (1996) Microbiol.
                     % doi:10.1099/13500872-142-8-1927 relatively stabile
                     % among growth rates, highest reported content is used
                     % to prevent overconstraining the model.
[f,~]       = measureAbundance(ecModel.enzymes);
GAM         = [];
pIDs        = data.genes;
sigma       = 1;

ecModel=setParam(ecModel,'lb','ATPM',0);


%% M145
sample      = [cellstr(repmat('M145',9,1));cellstr(repmat('M1152',8,1))];
sample(:,2) = num2cell([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]);
for i=1:9%length(sample(:,1))
    protData = data.mean(:,i) + data.std(:,i);
    prot(i).model = constrainEnzymes(ecModel,Ptot,sigma,f,[],pIDs,protData);%,GlcUptake);
    cd ../../..
    [prot(i).model, prot(i).gRate] = fixFluxes(prot(i).model,sample{i,1},sample{i,2});
    if prot(i).gRate<0
        prot(i).gRate=0.0070;
    end
    cd gecko/geckomat/limit_proteins
    [prot(i).model,prot(i).enzUsages,prot(i).modifications] = flexibilizeProteins(prot(i).model,prot(i).gRate);
end
