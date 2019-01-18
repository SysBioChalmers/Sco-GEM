%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = construct_ecScoGEM(model)
%
% Eduard Kerkhoven, 2018-12-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ecModel,model_data,kcats] = construct_ecScoGEM(model)

%Retrieve kcats & MWs for each rxn in model:
cd ../gecko/geckomat/get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,'streptomyces coelicolor');

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);

save('../../../../../scrap/ecScoGEM_noManualMods.mat','ecModel','model_data','kcats');
end