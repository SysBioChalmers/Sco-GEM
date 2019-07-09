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
save([root '/scrap/goodRxnsPool.mat'], 'goodRxnsPool*' );

cd('..')
for i=1:9
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},5000,true,false,false,goodRxnsPool1);
    [RS{i},rxns{i}] = organizeSolutions(model{i},RSraw{i},true,true);
end

for i=10:17
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},5000,true,false,false,goodRxnsPool2);
    [RS{i},rxns{i}] = organizeSolutions(model{i},RSraw{i},true,true);
end

for i=1:length(RSraw)
    tmpMean = mean(RS{i},2);
    tmpStd  = std(RS{i},0,2);
    [Lia, Locb] = ismember(rxns{1},rxns{i});
    RSmean(Lia,i) = tmpMean(Locb(Lia));
    RSsd(Lia,i) = tmpStd(Locb(Lia));
end

scatter(1:length(RSmean(:,1)),log10(RSmean(:,1)))
hold on
scatter(1:length(RSmean(:,10)),log10(RSmean(:,10)))
hold off

% Normalize for CO2 excretion
co2Idx = find(ismember(rxns{1},'EX_co2_e'));
for i=1:length(RSraw)
    normFactor = RSmean(co2Idx,i);
    disp(num2str(normFactor))
    RSmean(:,i) = RSmean(:,i) ./ normFactor;
    RSsd(:,i) = RSsd(:,i) ./ normFactor;
end

save([root '/scrap/RSpool.mat'],'RSraw','RSmean', 'RSsd','rxns')

%% Export data

% Calculate flux-Z
clear Z Zflux
Zflux(:,1) = getFluxZ(RSmean(:,2),RSmean(:,6));
Zflux(:,2) = getFluxZ(RSmean(:,2),RSmean(:,11));
Zflux(:,3) = getFluxZ(RSmean(:,11),RSmean(:,15));
Zflux(:,4) = getFluxZ(RSmean(:,6),RSmean(:,15));

rsOut=[rxns{1} num2cell(full(RSmean)) num2cell(full(RSsd)) num2cell(Zflux)];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RandSampComb_proteinPool_noProteomics.tsv'],'w');
fprintf(fid,[repmat('%s\t',1,38) '%s\n'],...
    ['rxns' strcat('MEAN_',transpose(sample)) strcat('STDEV_',transpose(sample)) 'Z_M145_29h_vs_M145_45h' 'Z_M145_29h_vs_M1152_41h' 'Z_M1152_41h_vs_M1152_57h' 'Z_M145_45h_vs_M1152_57h']);
for j=1:length(rxns{1})
    fprintf(fid,['%s\t' repmat('%d\t',1,37) '%d\n'],rsOut{j,:});
end
fclose(fid);
