%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ecModel = simulateProteomics(ecModel,data,sample)
%
% Eduard Kerkhoven, 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




load('gecko\models\scoGEM\EcscoGEM.mat')

[f,~] = measureAbundance(ecModel.enzymes);