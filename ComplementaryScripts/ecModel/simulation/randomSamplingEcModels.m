cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/proteomeModels.mat']);


%% Fix ATP maintenance before random sampling
for i=1:length(gRate)
    model{i} = setParam(model{i},'obj','ATPM',1);
    tmp = solveLP(model{i});
    disp(['Simulation ' num2str(i) ': ' num2str(tmp.f)])
    model{i} = setParam(model{i},'lb','ATPM',-0.95*tmp.f);
    model{i} = setParam(model{i},'ub','ATPM',-1.05*tmp.f);
end

%% First determine which reactions are part of loops and should therefore
% be discarded from further random sampling.
[~,goodRxns1] = randomSampling(model{1},2);
[~,goodRxns2] = randomSampling(model{11},2);
save([root '/scrap/goodRxns.mat'], 'goodRxns*' );

%% Take 1000 samples for each timepoint
for i=1:9
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},1000,true,false,false,goodRxns1);
    [RS{i},rxns] = organizeSolutions(model{i},RSraw{i});
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

for i=10:length(gRate)
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},1000,true,false,false,goodRxns2);
    [RS{i},rxns2] = organizeSolutions(model{i},RSraw{i});
    [Lia, Locb] = ismember(rxns2,rxns);
    tmpRS = zeros(length(rxns),1000);
    tmpRS(Locb(Lia),:) = sparse(RS{i}(Lia,:));
    RS{i} = tmpRS;
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

for i=1:length(RSraw)
    fid = fopen([root '/ComplementaryData/ecmodel/simulations/rsRaw_' sample{i} '.csv'],'w');
    for j=1:length(RSraw{i})
        fprintf(fid,'%s\t',model{i}.rxns{j});
        fprintf(fid,[repmat('%d\t',1,999) '%d\n'],full(RSraw{i}(j,:)));
    end
    fclose(fid);
end

RSmean=full(RSmean);
RSsd=full(RSsd);
save([root '/scrap/RScomb.mat'],'RS*','rxns')

for i=1:length(RS)
    [capUsage(:,i),absUsage(:,i)] = enzymeUsage(model{i},RSmean(:,i),true);
end

for i=1:length(RS)
    csvwrite([root '/ComplementaryData/ecmodel/simulations/RS_' num2str(i) '.csv'],full(RS{i}))
end

fid = fopen([root '/ComplementaryData/ecmodel/simulations/RS_rxns.csv'],'w');
fprintf(fid,'%s\n',string(rxns));
fclose(fid)
%% Write mean and standard deviations
rsOut=[rxns num2cell(RSmean) num2cell(RSsd)];
rsOut=[['' transpose(sample) strcat('STDEV_',transpose(sample))]; rsOut];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RandSamp.tab'],'w');
for k=1:length(rsOut)
   fprintf(fid,[repmat('%s\t',1,34) '%s\n'],rsOut(k,:));
end
fclose(fid);

%% Prepare for comparison with RNAseq/proteomics
rxngene=getRxnGeneList(modelScoGEM);
fid = fopen([root '/ComplementaryData/ecmodel/simulations/rxnGene.tab'],'w');
for k=1:length(rxngene)
   fprintf(fid,'%s\t%s\n',rxngene{k,:});
end
fclose(fid);

fid = fopen([root '/ComplementaryData/ecmodel/simulations/protein.tab'],'w');
tmp = ['gene\t' transpose(strcat(prot.sample,'\t'))];
fprintf(fid,[repmat('%s\t',1,51) '%s\n'],tmp{:});
for k=1:length(prot.norm)
   fprintf(fid,'%s\t',prot.genes{k});
   fprintf(fid,[repmat('%d\t',1,50) '%d\n'],prot.norm(k,:));
end
fclose(fid)