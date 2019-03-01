model=load ('https://github.com/SysBioChalmers/Sco-GEM/tree/fix/transporters/ModelFiles/xml/scoGEM_before_transporter_addition.mat');
Table1=readtable('Table1_updated_grRules_280119_300119_6.xlsx', 'ReadRowNames', false);
Table1=table2cell(Table1);
Table2=readtable('Table2_newTransportRxns_280119_300119_3.xlsx', 'ReadRowNames', false);
Table2=table2cell(Table2);
Table3=readtable('Table3_newTransportRxns_newMetabolites_280119_300119_3.xlsx', 'ReadRowNames', false);
Table3=table2cell(Table3);
Table4=readtable('ExRxns_added_Transporters_300119_3.xlsx', 'ReadRowNames', false);
Table4=table2cell(Table4);
Table5=readtable('Sink_Rxns_scoGEM_Transporters_300119_3.xlsx', 'ReadRowNames', false);
Table5=table2cell(Table5);
for i=1:1:27
    model=addReaction(model, Table2{i,1}, 'reactionFormula', Table2{i,3}, 'geneRule', Table2{i,2}, 'lowerBound', Table2{i,5}, 'upperBound',Table2{i,6});
end;
for i=1:1:23
    model=addReaction(model, Table3{i,1}, 'reactionFormula', Table3{i,3}, 'geneRule', Table3{i,2}, 'lowerBound', Table3{i,5}, 'upperBound', Table3{i,6});
end;
for i=1:1:23
    model=addReaction(model, Table4{i,1}, 'reactionFormula', Table4{i,2}, 'lowerBound', Table4{i,4}, 'upperBound', Table4{i,5});
end;
for i=1:1:23
    model=addReaction(model, Table5{i,1}, 'reactionFormula', Table5{i,2}, 'lowerBound', Table5{i,4}, 'upperBound', Table5{i,5});
end;
for i=1:1:37
   model=addReaction(model, Table1{i,1}, 'reactionFormula', Table1{i,3}, 'geneRule', Table1{i,2}, 'lowerBound', Table1{i,5}, 'upperBound',Table1{i,6});
end;
scoGEM_add_transporters=model;
scoGEM_add_transporters_change_dir_VALTA_LEUTA_ILETA = addReaction(scoGEM_add_transporters, 'VALTA', 'reactionFormula', 'akg[c] + val__L[c] -> glu__L[c] + 3mob[c]', 'lowerBound', 0, 'upperBound', 1000);
scoGEM_add_transporters_change_dir_VALTA_LEUTA_ILETA = addReaction(scoGEM_add_transporters_change_dir_VALTA_LEUTA_ILETA, 'LEUTA', 'reactionFormula', 'akg[c] + leu__L[c] -> glu__L[c] + 4mop[c] ', 'lowerBound', 0, 'upperBound', 1000);
scoGEM_model_added_transporters_cdir = addReaction(scoGEM_add_transporters_change_dir_VALTA_LEUTA_ILETA, 'ILETA', 'reactionFormula', 'akg[c] + ile__L[c] -> glu__L[c] + 3mop[c]', 'lowerBound', 0, 'upperBound', 1000);
scoGEM = scoGEM_model_added_transporters_cdir;
writeCbModel(scoGEM, 'fileName', 'scoGEM', 'format', 'sbml');


   

