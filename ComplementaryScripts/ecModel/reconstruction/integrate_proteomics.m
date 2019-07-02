%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrate_proteomics
%
% - Prepare proteomics data
% - Make sample-specific proteome constrained models
%
% Eduard Kerkhoven, 2019-03-27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/ecScoGEM_M145M1152.mat'])

%% - Prepare proteomics data
strain      = [string(repmat('M145',9,1));string(repmat('M1152',8,1))];
time        = string([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]);
sample      = strcat(strain,'_',time);

sigma   = 1;
Ptot    = 0.456; % As reported by Shahab et al. (1996) Microbiol.
cd([ecDir '/gecko/geckomat/limit_proteins'])
[f,~]   = measureAbundance(ecModel.enzymes);
cd(ecDir)
data    = prepareProteomicsM1152M145;
pIDs    = data.genes;
clear gRate modifications
%% - Make sample-specific proteome constrained models
cd([ecDir '/gecko/geckomat/limit_proteins'])
for i=1:length(sample)
    disp(['Generating ecModel for sample: ' sample(i)])
    protData = data.mean(:,i) + data.std(:,i);
    if i<10
        model{i} = constrainEnzymes(ecModel,Ptot,sigma,f,[],pIDs,protData);
    else
        model{i} = constrainEnzymes(ecModel_M1152,Ptot,sigma,f,[],pIDs,protData);
    end
    cd ([ecDir '\simulation'])
    [model{i}, gRate(i)] = fixFluxes(model{i},strain{i},time{i});
    if gRate(i)<0
        gRate(i)=0.0070;
    end
    cd([ecDir '/gecko/geckomat/limit_proteins'])
    [model{i},~,modifications{i}] = flexibilizeProteins(model{i},gRate(i));
    model{i} = setParam(model{i},'lb','BIOMASS_SCO_tRNA',0.95*gRate(i));
    model{i} = setParam(model{i},'ub','BIOMASS_SCO_tRNA',1.05*gRate(i));
    model{i} = setParam(model{i},'obj','ATPM',1);
end

%% All models should have the same flexibilization. First list
clear allModifications
allModifications = [modifications{1}.protein_IDs modifications{1}.modified_values];
for i=1:length(modifications)
    [Lia, Locb] = ismember(modifications{i}.protein_IDs,allModifications(:,1));
    c = {cell2mat(modifications{i}.modified_values(Lia));cell2mat(allModifications(Locb(Lia),2))};
    allModifications(Locb(Lia),2) = num2cell(max([c{:}], [], 2));
    allModifications = [allModifications;table2cell(modifications{i}(~Lia,[1,3]))];
end

for i=1:length(model)
    [Lia, Locb] = ismember(strcat(allModifications(:,1),'_exchange'),model{i}.rxnNames);
    model{i}.ub(Locb(Lia)) = cell2mat(allModifications(Lia,2));
end

for i=1:length(sample)
    model{i} = setParam(model{i},'obj','BIOMASS_SCO_tRNA',1);
    sol=solveLP(model{i},1);
    disp(['Sol = ' num2str(sol.f)])
end

save([root '/scrap/proteomeModels.mat'],'model','gRate','modifications','sample');

%% Export SBML files
for i=11:length(sample(:,1))
    exportModel(model{i},...
    [root,'/ModelFiles/xml/ec',sample{i},'.xml']);
end