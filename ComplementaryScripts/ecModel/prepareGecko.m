%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ecModel,model_data,kcats] = prepareGecko
%
% Eduard Kerkhoven, 20118-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepareGecko()

cd gecko/geckomat/get_enzyme_data

[swissprot,kegg] = updateDatabases;

save('../../databases/ProtDatabase.mat','swissprot','kegg');

end