cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/proteomeModels.mat']);

%% pFBA
clear fbaOut
cd(ecDir)
for i=1:length(gRate)
    modelTmp=model{i};
    modelTmp=setParam(modelTmp,'obj','ATPM',1);
    sol=solveLP(modelTmp,1);
    disp(['Simulation ' num2str(i) ': ' num2str(sol.f)])
    [solVect,rxns] = organizeSolutions(modelTmp,sol.x);
    fbaOut(:,i) = solVect;
end
fbaOut=[rxns num2cell(fbaOut)];
fbaOut=[['' transpose(sample)]; fbaOut];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-pFBA.tab'],'w');
for k=1:length(fbaOut)
   fprintf(fid,[repmat('%s\t',1,17) '%s\n'],fbaOut(k,:));
end
fclose(fid);
