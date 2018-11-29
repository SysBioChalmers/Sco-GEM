%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = constructEcScoGEM
%
% Eduard Kerkhoven, 20118-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = constructEcScoGEM(version)

%Set parameters, flags and options
format short e
name        = 'scoGEM';
org_name    = 'streptomyces coelicolor';
if nargin < 4
    version = datestr(now,'yymmdd'); % E.g. 181109 for 9 November 2018
end

%Load model
model       = importModel('../../ModelFiles/xml/scoGEM.xml');
sol=solveLP(model,1)

%Remove blocked rxns + correct model.rev
cd gecko/geckomat/change_model
[model,name,version] = preprocessModel(model,name,version);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);
save('../../../../../scrap/kcats.mat','model_data','kcats');

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

rates = dlmread('../../../../../ComplementaryData/growth/M145_estimated_rates.csv',';',1,0);

GlcUptake  = -rates(1,6);%
GluUptake  = -rates(1,4);%
gRate      = rates(1,2);%

ecModel=setParam(ecModel,'ub','EX_glc__D_e_REV',GlcUptake);
ecModel=setParam(ecModel,'ub','EX_glu__L_e_REV',GluUptake);
ecModel=setParam(ecModel,'ub','EX_nh4_e_REV',0);
%ecModel=setParam(ecModel,'lb','ATPM',2.64);


save(['../../../../../scrap/Ec' name '.mat'],'ecModel','modifications') 
load(['../../scrap/Ec' name '.mat']) 
%Constrain model to batch conditions:

sigma       = 0.40;
Ptot        = 0.429; %Assumed constant. Taken as average from growth rates
                     %between0 0.024 and 0.195, as reported by Shabab et
                     %al. (1996) Microbiol. doi:10.1099/13500872-142-8-1927 
cd ../limit_proteins

[ecModel_pool,OptSigma] = getConstrainedModel(ecModel,'',sigma,Ptot,gRate,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
save(['../../../../../scrap/Ec' name '_pool.mat'],'ecModel_pool','modifications') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%