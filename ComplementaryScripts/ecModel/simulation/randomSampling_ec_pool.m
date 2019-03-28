cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/ecPool.mat'],'ecModel*_pool')

strain      = [string(repmat('M145',9,1));string(repmat('M1152',8,1))];
time        = string([21;29;33;37;41;45;49;53;57;33;41;45;49;53;57;61;65]);
sample      = strcat(strain,'_',time);

f = 0.4247;
sigma   = 0.5;
Ptot    = 0.456;
ecModel_M145_pool = setParam(ecModel_M145_pool,'ub','prot_pool_exchange',f*sigma*Ptot);
ecModel_M1152_pool = setParam(ecModel_M1152_pool,'ub','prot_pool_exchange',f*sigma*Ptot);

for i=1:length(sample)
    disp(['Generating ecModel for sample: ' sample{i}])
    cd ([ecDir '\simulation'])
    if i<10
        [model{i}, gRate(i)] = fixFluxes(ecModel_M145_pool,strain{i},time{i});
    else
        [model{i}, gRate(i)] = fixFluxes(ecModel_M1152_pool,strain{i},time{i});
    end
    if gRate(i)<0
        gRate(i)=0.0070;
    end
end

%% Random sampling
for i=1:length(gRate)
    model{i} = setParam(model{i},'lb','BIOMASS_SCO_tRNA',0.99*gRate(i));
    model{i} = setParam(model{i},'ub','BIOMASS_SCO_tRNA',1.01*gRate(i));
    model{i} = setParam(model{i},'obj','ATPM',1);
    tmp = solveLP(model{i});
    disp(['ATP maintenance reached in ' num2str(i) ': ' num2str(-tmp.f)])
    model{i} = setParam(model{i},'lb','ATPM',-0.99*tmp.f);
    %model{i} = setParam(model{i},'ub','ATPM',-1.01*tmp.f);
    model{i} = setParam(model{i},'obj','prot_pool_exchange',-1);
    tmp = solveLP(model{i});
    disp(['Minimum protein required in ' num2str(i) ': ' num2str(-tmp.f)])
    model{i} = setParam(model{i},'lb','prot_pool_exchange',0.99*tmp.f);
    model{i} = setParam(model{i},'ub','prot_pool_exchange',1.01*tmp.f);
end

[~,goodRxnsPool1] = randomSampling(model{1},2);
[~,goodRxnsPool2] = randomSampling(model{11},2);
load([root '/scrap/goodRxnsPool.mat'], 'goodRxnsPool*' );

cd('..')
for i=1:9
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},1000,true,false,false,goodRxnsPool1);
    [RS{i},rxns] = organizeSolutions(model{i},RSraw{i},true,true);
end

for i=10:17
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},1000,true,false,false,goodRxnsPool2);
    [RS{i},rxns2] = organizeSolutions(model{i},RSraw{i},true,true);
    [Lia, Locb] = ismember(rxns2,rxns);
    tmpRS = zeros(length(rxns),1000);
    tmpRS(Locb(Lia),:) = sparse(RSraw{i}(Lia,:));
    RS{i} = tmpRS;
end

% Normalize for CO2 excretion
co2Idx = find(ismember(rxnsRaw,'EX_co2_e'));
for i=1:length(RS)
    normFactor = RS{i}(co2Idx,:);
    disp(num2str(normFactor(1)))
    RSnorm{i} = RS{i} ./ normFactor;
    RSmean(:,i) = mean(RSnorm{i},2);
    RSsd(:,i) = std(RSnorm{i},0,2);
end

%% Export data
save([root '/scrap/RSpool.mat'],'RSraw','RS','rxns')

for i=1:length(RS)
    RSmean(:,i) = mean(RSnorm{i},2);
    RSsd(:,i) = std(RSnorm{i},0,2);
end

% Calculate flux-Z
clear Z Zflux
Zflux(:,1) = getFluxZ(RSnorm{2},RSnorm{6});
Zflux(:,2) = getFluxZ(RSnorm{2},RSnorm{11});
Zflux(:,3) = getFluxZ(RSnorm{11},RSnorm{15});
Zflux(:,4) = getFluxZ(RSnorm{6},RSnorm{15});

rsOut=[rxns num2cell(RSmean) num2cell(RSsd) num2cell(Zflux)];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RSPoolCombNorm.tsv'],'w');
fprintf(fid,[repmat('%s\t',1,38) '%s\n'],...
    ['rxns' strcat('MEAN_',transpose(sample)) strcat('STDEV_',transpose(sample)) 'Z_M145_29h_vs_M145_45h' 'Z_M145_29h_vs_M1152_41h' 'Z_M1152_41h_vs_M1152_57h' 'Z_M145_45h_vs_M1152_57h']);
for j=1:length(rxnsRaw)

    fprintf(fid,['%s\t' repmat('%d\t',1,37) '%d\n'],rsOut{j,:});
end
fclose(fid);

% RxnGene List
modelScoGEM = importModel('../../ModelFiles/xml/scoGEM.xml');
rxngene=getRxnGeneList(modelScoGEM);
fid = fopen([root '/ComplementaryData/ecmodel/simulations/rxnGenePool.tsv'],'w');
for k=1:length(rxngene)
    fprintf(fid,'%s\t%s\t',rxngene{k,:});
    fprintf(fid,'%s\n',modelScoGEM.rxnNames{ismember(modelScoGEM.rxns,rxngene{k,1})});
end
fclose(fid);
