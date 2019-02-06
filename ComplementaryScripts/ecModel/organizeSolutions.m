function [solVect,rxns] = organizeSolutions(model,sol,verbose);

if nargin<3
    verbose = false;
end
%% Combine isoenzymatic reactions
rxns = model.rxns;
keep = ones(length(rxns),1);
for i=1:length(rxns)
    if keep(i) == 1
        if contains(rxns(i),'arm_')
            rxnId = regexprep(rxns{i},'arm_','');
            rxns{i} = rxnId;
            keep(find(contains(rxns,[rxnId 'No']))) = 0;
        end
        rxns{i} = regexprep(rxns{i},'No1$','');
        if contains(rxns(i),'_REV')
            rxnId = regexprep(rxns{i},'_REV','');
            fwdId = find(ismember(rxns,rxnId));
            if isempty(fwdId)
                if verbose
                    disp(['No fwd reaction for ' rxnId ' found'])
                end
                sol(i) = -sol(i);
            else
                sol(fwdId) = sol(fwdId) - sol(i);
            end
            keep(i) = 0;
        end
    end
end
keep=logical(keep);
rxns = rxns(keep);
sol  = sol(keep);

%% Organize enzyme utilization reactions
%Make list of enzymes in the model
proteins  = sort(model.enzymes);
%Index of first enzyme
protStart = find(contains(rxns,'prot_'));
protStart = protStart(1);
%Trim down to only protein identifier
rxnNames  = rxns(protStart:end-1);
rxnNames  = regexprep(rxnNames,'(draw_)?prot_','');
rxnNames  = regexprep(rxnNames,'_exchange','');
%Match enzyme reactions to protein list
[~,idx]=ismember(proteins,rxnNames);
idx = idx + protStart - 1;
%Reorder solution vector
sol(protStart:end-1) = sol(idx);
solVect = sol;
%Write vector specifying order of reactions in solution vector
%rxns(1:protStart-1,1) = model.rxns(1:protStart-1);
proteins=strcat('prot_',proteins);
rxns(protStart:end-1) = proteins;
%rxns(end+1) = model.rxns(end);

end

