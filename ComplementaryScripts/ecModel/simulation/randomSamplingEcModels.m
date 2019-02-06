cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/proteomeModels.mat']);

%% Fix ATP maintenance before random sampling
for i=1:length(gRate)
    tmp = solveLP(model{i});
    disp(['Simulation ' num2str(i) ': ' num2str(tmp.f)])
    model{i} = setParam(model{i},'lb','ATPM',-0.99*tmp.f);
    model{i} = setParam(model{i},'ub','ATPM',-1.01*tmp.f);
end

%% First determine which reactions are part of loops and should therefore
% be discarded from further random sampling.
[~,goodRxns1] = randomSampling(model{1},2);
[~,goodRxns2] = randomSampling(model{11},2);
save([root '/scrap/goodRxns.mat'], 'goodRxns*' );

%% Take 1000 samples for each timepoint
for i=1:9
    disp(['Sample: ' num2str(i)])
    RS{i} = randomSampling(model{i},1000,true,false,false,goodRxns1);
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

for i=10:length(gRate)
    disp(['Sample: ' num2str(i)])
    RS{i} = randomSampling(model{i},1000,true,false,false,goodRxns2);
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

RSmean=full(RSmean);
RSsd=full(RSsd);
save([root '/scrap/RS.mat'],'RS*')

%% Write mean and standard deviations
clear mnOut sdOut
for i=1:length(gRate)
    [mnVect,rxns] = organizeSolutions(model{i},RSmean(:,i));
    [sdVect,rxns] = organizeSolutions(model{i},RSsd(:,i));
    mnOut(:,i) = mnVect;
    sdOut(:,i) = sdVect;
end

rsOut=[rxns num2cell(mnOut) num2cell(sdOut)];
rsOut=[['' transpose(sample) strcat('STDEV_',transpose(sample))]; rsOut];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RandSamp.tab'],'w');
for k=1:length(rsOut)
   fprintf(fid,[repmat('%s\t',1,34) '%s\n'],rsOut(k,:));
end
fclose(fid);