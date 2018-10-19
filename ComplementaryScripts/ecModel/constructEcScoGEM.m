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
    version = '0.1';
end

%Load model
model       = importModel('../ModelFiles/xml/scoGEM.xml');

%Remove blocked rxns + correct model.rev and reorder biomass:
cd ../gecko/change_model
[model,name,version] = preprocessModel(model,name,version);
model                = reorderBiomass(model);

%Set reaction to irreversible
model=setParam(model,'lb','BFBP',0); % fructose bisphophatase

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);
save(['../../models/' name '/data/' name '_enzData_allUniprot.mat'],'model_data','kcats','version')

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
[ecModel,modifications] = manualModifications(ecModel);
% Modify biomass
cd ../../ComplementaryData/curation

ecModel=setParam(ecModel,'eq','EX_glc__D_e',0);
ecModel=setParam(ecModel,'ub','EX_glc__D_e_REV',0.51);
ecModel=setParam(ecModel,'obj','BM_growth',1);
ecModel=setParam(ecModel,'lb','ATPM',3);
sol=solveLP(ecModel)


%Constrain model to batch conditions:
sigma    = 0.4;      %Optimized for glucose
Ptot     = 0.4120;    %Assumed constant
gR_exp=0.03427;
c_source = 'D-glucose exchange (reversible)'; %Rxn name for the glucose uptake reaction
cd ../limit_proteins
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,gR_exp,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])

%Save output models:
cd ../../models
ecModel = saveECmodel(ecModel,toolbox,name,version);
saveECmodel(ecModel_batch,toolbox,[name '_batch'],version);
cd ../geckomat

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%