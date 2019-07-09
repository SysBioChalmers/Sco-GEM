cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/proteomeModels.mat']);

%% Fix ATP maintenance before random sampling
for i=1:length(gRate)
    model{i} = setParam(model{i},'obj','ATPM',1);
    tmp = solveLP(model{i});
    disp(['Simulation ' num2str(i) ': ' num2str(tmp.f)])
    model{i} = setParam(model{i},'lb','ATPM',-0.99*tmp.f);
    model{i} = setParam(model{i},'ub','ATPM',-1.01*tmp.f);
end

%% First determine which reactions are part of loops and should therefore
% be discarded from further random sampling.
[~,goodRxns1] = randomSampling(model{1},2);
[~,goodRxns2] = randomSampling(model{11},2);
save([root '/scrap/goodRxns2.mat'], 'goodRxns*' );

%% Take 5000 samples for each timepoint. Reorder reactions so that they
% match between models.
nSamples = 5000;
for i=1:9
    disp(['Sample: ' num2str(i)])
    %RSraw{i} = randomSampling(model{i},nSamples,true,false,false,goodRxns1);
    [RS{i},rxns{i}] = organizeSolutions(model{i},RSraw{i},true,true);
end

for i=10:17
    disp(['Sample: ' num2str(i)])
    %RSraw{i} = randomSampling(model{i},nSamples,true,false,false,goodRxns2);
    [RS{i},rxns{i}] = organizeSolutions(model{i},RSraw{i},true,true);
end

% Export data
save([root '/scrap/RScomb.mat'],'RSraw','RS','rxns')

for i=1:length(RSraw)
    tmpMean = mean(RS{i},2);
    tmpStd  = std(RS{i},0,2);
    [Lia, Locb] = ismember(rxns{1},rxns{i});
    RSmean(Lia,i) = tmpMean(Locb(Lia));
    RSsd(Lia,i) = tmpStd(Locb(Lia));
end

% Normalize for CO2 excretion
co2Idx = find(ismember(rxns{1},'EX_co2_e'));
for i=1:length(RSraw)
    normFactor = RSmean(co2Idx,i);
    disp(num2str(normFactor))
    RSmean(:,i) = RSmean(:,i) ./ normFactor;
    RSsd(:,i) = RSsd(:,i) ./ normFactor;
end

clear Z Zflux
Zflux(:,1) = getFluxZ(RSmean(:,2),RSmean(:,6));
Zflux(:,2) = getFluxZ(RSmean(:,2),RSmean(:,11));
Zflux(:,3) = getFluxZ(RSmean(:,11),RSmean(:,15));
Zflux(:,4) = getFluxZ(RSmean(:,6),RSmean(:,15));

rsOut=[rxns{1} num2cell(full(RSmean)) num2cell(full(RSsd)) num2cell(Zflux)];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RandSampComb_CO2norm.tsv'],'w');
fprintf(fid,[repmat('%s\t',1,39) '%s\n'],...
    ['rxns' strcat('MEAN_',transpose(sample)) strcat('STDEV_',transpose(sample)) 'Z_M145_29h_vs_M145_45h' 'Z_M145_29h_vs_M1152_41h' 'Z_M1152_41h_vs_M1152_57h' 'Z_M145_45h_vs_M1152_57h']);
for j=1:length(rsOut)
    fprintf(fid,['%s' repmat('\t%d',1,38) '\n'],rsOut{j,:});
end
fclose(fid);