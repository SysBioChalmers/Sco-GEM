%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = constructEcScoGEM
%
% Eduard Kerkhoven, 2018-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = constructEcScoGEM

%Load model
model       = importModel('../../ModelFiles/xml/scoGEM.xml');
sol=solveLP(model,1);
fprintf(['Growth rate of template model: ' num2str(-sol.f) '\n']);

%Remove blocked rxns + correct model.rev
cd reconstruction
model = preprocessModel(model);

%Retrieve kcats & MWs for each rxn in model:
cd ../gecko/geckomat/get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,'streptomyces coelicolor');
save('../../../../../ComplementaryData/ecmodel/kcats.mat','model_data','kcats');
load('../../../../../ComplementaryData/ecmodel/kcats.mat');

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
cd ../../../reconstruction
[ecModel,modifications] = manualModifications(ecModel);

rates = dlmread('../../../ComplementaryData/growth/M145_estimated_rates.csv',';',1,0);

GlcUptake  = -rates(1,4);%
GluUptake  = -rates(1,5);%
gRate      = rates(1,3);%

ecModel=setParam(ecModel,'ub','EX_glc__D_e_REV',GlcUptake);
ecModel=setParam(ecModel,'ub','EX_glu__L_e_REV',GluUptake);
ecModel=setParam(ecModel,'ub','EX_nh4_e_REV',0);

save('../../../ModelFiles/mat/ecScoGEM.mat','ecModel','modifications')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%