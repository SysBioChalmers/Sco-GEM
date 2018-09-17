%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = addSBOterms(model)
%
%   Adds relevant SBO terms to the model. Uses RAVEN 2.0.3+
%
% 2018-09-15    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = addSBOterms(model)

% Add SBO terms for mets. No 'pseudometabolites' exist, so all are
% annotated with SBO:00000247, 'simple chemical'. Pseudometabolite should
% likely be annotated with SB:0000649, 'biomass'.
model=addMiriams(model,'mets',model.mets,'SBO','0000247');

% Add SBO terms for rxns
for i=1:length(model.rxns)
    rxnName   = model.rxnNames{i};
    metNames  = model.metNames(model.S(:,i) ~= 0);
    metComps  = model.metComps(model.S(:,i) ~= 0);
    metStoich = model.S(model.S(:,i) ~= 0,i);
    if contains(rxnName,'biomass','IgnoreCase',true)
        % Biomass production
        model=addMiriams(model,'rxns',model.rxns(i),'SBO','0000629');
    elseif contains(rxnName,'maintenance','IgnoreCase',true)
        % ATP maintenance
        model=addMiriams(model,'rxns',model.rxns(i),'SBO','0000630');
    elseif numel(metNames)~=numel(unique(metNames))>0
        % Transport reaction
        model=addMiriams(model,'rxns',model.rxns(i),'SBO','0000655');
    elseif length(metNames) == 1
        if strcmp(model.comps{metComps},'e')
            % Exchange reaction
            model=addMiriams(model,'rxns',model.rxns(i),'SBO','0000627');
        else
            % Sink reaction
            model=addMiriams(model,'rxns',model.rxns(i),'SBO','0000632');
        end
    end
end
end
