%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = constructEcScoGEM
%
% Eduard Kerkhoven, 2018-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = constructEcScoGEM(version)

%Set parameters, flags and options
%format short e
name        = 'scoGEM';
org_name    = 'streptomyces coelicolor';
if nargin < 4
    version = datestr(now,'yymmdd'); % E.g. 181109 for 9 November 2018
end

%Load model
model       = importModel('../../ModelFiles/xml/scoGEM.xml');
sol=solveLP(model,1);
fprintf(['Growth rate of template model: ' num2str(-sol.f) '\n']);

%Remove blocked rxns + correct model.rev
cd gecko/geckomat/change_model
[model,name,version] = preprocessModel(model,name,version);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);
save('../../../../../ComplementaryData/ecmodel/kcats.mat','model_data','kcats');

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

rates = dlmread('../../../../../ComplementaryData/growth/M145_estimated_rates.csv',';',1,0);

GlcUptake  = -rates(1,4);%
GluUptake  = -rates(1,5);%
gRate      = rates(1,3);%

ecModel=setParam(ecModel,'ub','EX_glc__D_e_REV',GlcUptake);
ecModel=setParam(ecModel,'ub','EX_glu__L_e_REV',GluUptake);
ecModel=setParam(ecModel,'ub','EX_nh4_e_REV',0);
%ecModel=setParam(ecModel,'lb','ATPM',2.64);

save(['../../../../../ModelFiles/mat/Ec' name '.mat'],'ecModel','modifications')
load(['../../../../../ModelFiles/mat/Ec' name '.mat'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%