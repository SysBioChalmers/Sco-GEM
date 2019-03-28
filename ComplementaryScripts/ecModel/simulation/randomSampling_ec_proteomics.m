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
load([root '/scrap/goodRxns2.mat'], 'goodRxns*' );

%% Take 1000 samples for each timepoint. Reorder reactions so that they
% match between models.
for i=2:9
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},1000,true,false,false,goodRxns1);
    cd([ecDir '/gecko/geckomat/utilities'])
    RSmean(:,i) = mean(RSraw{i},2);
    [capUsage{i},absUsage{i}] = enzymeUsage(model{i},RSmean(:,i));
    cd(ecDir)
    [RS{i},rxnsRaw] = organizeSolutions(model{i},RSraw{i},false,true);
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

for i=10:17
    disp(['Sample: ' num2str(i)])
    RSraw{i} = randomSampling(model{i},1000,true,false,false,goodRxns2);
    cd([ecDir '/gecko/geckomat/utilities'])
    tmpMean = mean(RSraw{i},2);
    [capUsage{i},absUsage{i}] = enzymeUsage(model{i},tmpMean);
    cd(ecDir)
    [Lia, Locb] = ismember(model{i}.rxns,rxnsRaw);
    tmpRS = zeros(length(rxnsRaw{1}),1000);
    tmpRS(Locb(Lia),:) = sparse(RSraw{i}(Lia,:));
    RS{i} = tmpRS;
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

% Export data
RSmeanRaw=full(RSmean);
RSsdRaw=full(RSsd);
load([root '/scrap/RScomb.mat'],'RS*','rxnsRaw')

for i=1:length(RSraw)
    fid = fopen([root '/ComplementaryData/ecmodel/simulations/rsRaw_' sample{i} '.tsv'],'w');
    for j=1:length(RSraw{i})
        fprintf(fid,'%s\t',rxnsRaw{j});
        fprintf(fid,[repmat('%d\t',1,999) '%d\n'],full(RSraw{i}(j,:)));
    end
    fclose(fid);
end

% Calculate flux-Z
clear Z Zflux
Zflux(:,1) = getFluxZ(RSraw{2},RSraw{6});
Zflux(:,2) = getFluxZ(RSraw{2},RSraw{11});
Zflux(:,3) = getFluxZ(RSraw{11},RSraw{15});
Zflux(:,4) = getFluxZ(RSraw{6},RSraw{15});
fid = fopen([root '/ComplementaryData/ecmodel/simulations/fluxZ.tsv'],'w');
for k=1:length(Zflux)
   fprintf(fid,'%s\t',rxnsRaw{k});
   fprintf(fid,'%d\t%d\t%d\t%d\n',Zflux(k,:));
end
fclose(fid);


rsOut=[rxnsRaw num2cell(RSmeanRaw) num2cell(RSsdRaw) num2cell(Zflux)];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RandSamp.tsv'],'w');
fprintf(fid,[repmat('%s\t',1,38) '%s\n'],...
    ['rxns' strcat('MEAN_',transpose(sample)) strcat('STDEV_',transpose(sample)) 'Z_M145_29h_vs_M145_45h' 'Z_M145_29h_vs_M1152_41h' 'Z_M1152_41h_vs_M1152_57h' 'Z_M145_45h_vs_M1152_57h']);
for j=1:length(rxnsRaw)
    fprintf(fid,['%s\t' repmat('%d\t',1,37) '%d\n'],rsOut{j,:});
end
fclose(fid);


% RxnGene List
rxngene=getRxnGeneList(model{1});
fid = fopen([root '/ComplementaryData/ecmodel/simulations/rxnGene.tsv'],'w');
for k=1:length(rxngene)
   fprintf(fid,'%s\t%s\n',rxngene{k,:});
end
fclose(fid);

% Proteomics data
prot=prepareProteomicsM1152M145();
rsOut=[prot.genes num2cell(prot.norm)];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/proteomeNorm.tsv'],'w');
fprintf(fid,[repmat('%s\t',1,51) '%s\n'],['genes' string(transpose(prot.sample))]);
for j=1:length(prot.genes)
    fprintf(fid,['%s\t' repmat('%d\t',1,50) '%d\n'],rsOut{j,:});
end
fclose(fid);

for j=1:length(rsOut)
    fprintf(fid,['%s\t' repmat('%d\t',1,33) '%d\n'],rsOut{j,:});
end
fclose(fid);

% Protein Z
clear Z
Z(:,1) = getFluxZ(prot.norm(:,[2,11,20]),prot.norm(:,[6,15,24]));
Z(:,2) = getFluxZ(prot.norm(:,[2,11,20]),prot.norm(:,[29,37,45]));
Z(:,3) = getFluxZ(prot.norm(:,[29,37,45]),prot.norm(:,[33,41,49]));
Z(:,4) = getFluxZ(prot.norm(:,[6,15,24]),prot.norm(:,[33,41,49]));
fid = fopen([root '/ComplementaryData/ecmodel/simulations/protZ.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["genes","EL145","EE","EL1152","LL"]);
for k=1:length(Z)
   fprintf(fid,'%s\t',prot.genes{k});
   fprintf(fid,'%d\t%d\t%d\t%d\n',Z(k,:));
end
fclose(fid);

%% Calculate for contracted model
clear RSmean RSsd
for i=1:length(RSraw)
	disp(['Sample: ' num2str(i)])
    if i < 10
        [RS{i},rxns] = organizeSolutions(model{i},RSraw{i},true);
    else
        [RS{i},rxns2] = organizeSolutions(model{i},RSraw{i},true);
        [Lia, Locb] = ismember(rxns2,rxns);
        rxns2(Lia);
        tmpRS = zeros(length(rxns),100);
        tmpRS(Locb(Lia),:) = sparse(RS{i}(Lia,:));
        RS{i} = tmpRS;
    end
    RSmean(:,i) = mean(RS{i},2);
    RSsd(:,i) = std(RS{i},0,2);
end

co2Idx = find(ismember(rxns,'EX_co2_e'));
for i=1:length(RS)
    normFactor = RS{i}(co2Idx,:);
    disp(num2str(normFactor))
    RSNorm{i} = RS{i} ./ normFactor;
    RSmean(:,i) = mean(RSNorm{i},2);
    RSsd(:,i) = std(RSNorm{i},0,2);
end

clear Z Zflux
Zflux(:,1) = getFluxZ(RS{2},RS{6});
Zflux(:,2) = getFluxZ(RS{2},RS{11});
Zflux(:,3) = getFluxZ(RS{11},RS{15});
Zflux(:,4) = getFluxZ(RS{6},RS{15});
fid = fopen([root '/ComplementaryData/ecmodel/simulations/fluxZnorm.tsv'],'w');
for k=1:length(Zflux)
   fprintf(fid,'%s\t',rxns{k});
   fprintf(fid,'%d\t%d\t%d\t%d\n',Zflux(k,:));
end
fclose(fid);

modelScoGEM = importModel('../../ModelFiles/xml/scoGEM.xml');
[a,b] = ismember(rxns,modelScoGEM.rxns);

emp = cell(numel(rxns),1);
emp(:) = {''};
grRules = emp;
grRules(a) = modelScoGEM.grRules(b(a));

rxnNames = emp;
rxnNames(a) = modelScoGEM.rxnNames(b(a));

rsOut=[rxns rxnNames grRules num2cell(full(RSmean)) num2cell(full(RSsd)) num2cell(Zflux)];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-RandSampComb.tsv'],'w');
fprintf(fid,[repmat('%s\t',1,40) '%s\n'],...
    ['rxns' 'rxnNames' 'grRules' strcat('MEAN_',transpose(sample)) strcat('STDEV_',transpose(sample)) 'Z_M145_29h_vs_M145_45h' 'Z_M145_29h_vs_M1152_41h' 'Z_M1152_41h_vs_M1152_57h' 'Z_M145_45h_vs_M1152_57h']);
for j=1:length(rsOut)
    fprintf(fid,['%s\t%s\t%s\t' repmat('%d\t',1,37) '%d\n'],rsOut{j,:});
end
fclose(fid);

for i=1:length(RS)
    fid = fopen([root '/ComplementaryData/ecmodel/simulations/rsComb_' sample{i} '.tsv'],'w');
    for j=1:length(RS{i})
        fprintf(fid,'%s\t',rxns{j});
        fprintf(fid,[repmat('%d\t',1,999) '%d\n'],full(RS{i}(j,:)));
    end
    fclose(fid);
end




%% Prepare for comparison with RNAseq/proteomics
rxngene=getRxnGeneList(model{1});
fid = fopen([root '/ComplementaryData/ecmodel/simulations/rxnGene.tab'],'w');
for k=1:length(rxngene)
   fprintf(fid,'%s\t%s\n',rxngene{k,:});
end
fclose(fid);
fid = fopen([root '/ComplementaryData/ecmodel/simulations/rxns.tab'],'w');
for k=1:length(model{1}.rxns,rxngene)
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