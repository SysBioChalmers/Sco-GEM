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

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);

ecModel=setParam(ecModel,'eq','EX_glc__D_e',0);
ecModel=setParam(ecModel,'ub','EX_glc__D_e_REV',0.51);
ecModel=setParam(ecModel,'lb','ATPM',0);
sol=solveLP(ecModel)

save(['../../models/Ec' name '.mat'],'ecModel','modifications') 
%Constrain model to batch conditions:
sigma    = 0.6;      %Optimized for glucose
Ptot     = 0.4120;    %Assumed constant
obj_Val=0.03427;
c_source = 'D-glucose exchange (reversible)'; %Rxn name for the glucose uptake reaction
cd ../limit_proteins

[ecModel_pool,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,obj_Val,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

save('gecko/models/ecScoGEM_pool.mat','ecModel_pool')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%