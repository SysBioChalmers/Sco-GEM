cd ..
ecDir = pwd();
root  = regexprep(ecDir,'(.*)\\[^\\]*\\.*','$1');
load([root '/scrap/proteomeModels.mat']);

%% pFBA
clear fbaOut
for i=1:length(gRate)
    modelTmp=model{i};
    sol=solveLP(modelTmp,1);
    [solVect,rxns] = organizeSolutions(modelTmp,sol.x);
    fbaOut(:,i) = solVect;
    disp(['Simulation ' num2str(i) ': ' num2str(sol.f)])
end
fbaOut=[rxns num2cell(fbaOut)];
fbaOut=[['' transpose(sample)]; fbaOut];
fid = fopen([root '/ComplementaryData/ecmodel/simulations/ec-pFBA.tab'],'w');
for k=1:length(fbaOut)
   fprintf(fid,[repmat('%s\t',1,17) '%s\n'],fbaOut(k,:));
end
fclose(fid);
