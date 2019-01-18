%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
%
% Benjamin J. Sanchez. Last edited: 2017-10-29
% Ivan Domenzain.      Last edited: 2018-05-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,modifications] = manualModifications(model)

%Read manual data:
fID           = fopen('../gecko/databases/manual_data.txt');
data          = textscan(fID,'%s %s %s %s %f','delimiter','\t');
structure     = data{2};
protGenes     = data{4};
kcats         = data{5}.*3600;
data          = load('../gecko/databases/ProtDatabase.mat');
swissprot     = data.swissprot;
kegg          = data.kegg;
fclose(fID);
modifications{1} = cell(0,1);
modifications{2} = cell(0,1);

%Construct curated complexes:
uniprots = cell(size(kcats));
stoich   = cell(size(kcats));
for i = 1:length(kcats)
    uniprots{i}  = strsplit(structure{i},' + ');
    stoich{i}    = ones(size(uniprots{i}));
    %Separate complex strings in units and amount of each unit:
    for j = 1:length(uniprots{i})
        unit = uniprots{i}{j};
        pos  = strfind(unit,' ');
        if isempty(pos)
            stoich{i}(j)   = 1;
            uniprots{i}{j} = unit;
        else
            stoich{i}(j)   = str2double(unit(1:pos-1));
            uniprots{i}{j} = unit(pos+1:end);
        end
    end
end

for i = 1:length(model.rxns)
    reaction = model.rxnNames{i};
    %Find set of proteins present in rxn:
    S        = model.S;
    subs_pos = find(S(:,i) < 0);
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos);
    prot_set = cell(size(int_pos));
    MW_set   = 0;
    for j = 1:length(int_pos)
        met_name    = model.mets{int_pos(j)};
        prot_set{j} = met_name(6:end);
        MW_set      = MW_set + model.MWs(strcmp(model.enzymes,prot_set{j}));
    end
    %Find intersection with manual curated data:
    for j = 1:length(uniprots)
        int    = intersect(prot_set,uniprots{j});
        if length(int)/max(length(prot_set),length(uniprots{j})) > 0.50 % 50% match
            %Erase previous protein stoich. coeffs from rxn:
            for k = 1:length(prot_set)
                model.S(int_pos(k),i) = 0;
            end
            %If some proteins where not present previously, add them:
            newMets = uniprots{j};
            grRule  = protGenes{j};
            for k = 1:length(uniprots{j})
                if sum(strcmp(model.enzymes,uniprots{j}{k})) == 0
                    model = addProtein(model,uniprots{j}{k},kegg,swissprot);
                end
                newMets{k} = ['prot_' newMets{k}];
            end
            %Add new protein stoich. coeffs to rxn:
            kvalues = kcats(j)./stoich{j};
            rxnID   = model.rxns{i};
            rxnName = model.rxnNames{i};
            model   = addEnzymesToRxn(model,kvalues,rxnID,newMets,{rxnID,rxnName},grRule);
        end
    end
    %Update int_pos:
    S        = model.S;
    subs_pos = find(S(:,i) < 0);
    %Get the proteins that are part of the i-th rxn
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = transpose(intersect(subs_pos,prot_pos));
    %%%%%%%%%%%%%%%%%%%%%%%%%%  Individual Changes:  %%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(int_pos)
        enzName = model.mets(int_pos(j));
        %%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:
        [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications);
        if ~isempty(newValue)
            model.S(int_pos(j),i) = newValue;
        end
    end
    if rem(i,100) == 0 | i == rxns
        disp(['Improving model with curated data: Ready with rxn ' num2str(i)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove repeated reactions (2017-01-16):
rem_rxn = false(size(model.rxns));
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) && model.lb(i) == model.lb(j) && ...
                model.ub(i) == model.ub(j)
            rem_rxn(j) = true;
            disp(['Removing repeated rxn: ' model.rxns{i} ' & ' model.rxns{j}])
        end
    end
end
model = removeReactions(model,model.rxns(rem_rxn),true,true);
% Merge arm reactions to reactions with only one isozyme (2017-01-17):
arm_pos = zeros(size(model.rxns));
p       = 0;
for i = 1:length(model.rxns)
    rxn_id = model.rxns{i};
    if contains(rxn_id,'arm_')
        rxn_code = rxn_id(5:end);
        k        = 0;
        for j = 1:length(model.rxns)
            if ~isempty(strfind(model.rxns{j},[rxn_code 'No']))
                k      = k + 1;
                pos    = j;
                grRule = model.grRules{j};
            end
        end
        if k == 1
            %Condense both reactions in one:
            equations.mets         = model.mets;
            equations.stoichCoeffs = model.S(:,i) + model.S(:,pos);
            model = changeRxns(model,model.rxns(pos),equations);
            model.grRules{pos} = grRule;
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end

% Aconitase (Q7AKF3/EC 4.2.1.3): The rxn is represented as a two
% step rxn, so the kcat must be divided by 2
index = find(strcmpi('prot_Q7AKF3',model.mets));
rxnIndxs = find(model.S(index,:));
rxnIndxs = rxnIndxs(1:end-2);
model.S(index,rxnIndxs) = model.S(index,rxnIndxs)/2;

% Remove saved arm reactions:
model = removeReactions(model,model.rxns(arm_pos(1:p)),true,true);

% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_enz(i) = true;
    end
end
rem_enz = model.enzymes(rem_enz);
for i = 1:length(rem_enz)
    model = deleteProtein(model,rem_enz{i});
    disp(['Removing unused protein: ' rem_enz{i}])
end

% Remove incorrect pathways:
cd ../gecko/geckomat/change_model/
model = removeIncorrectPathways(model);
cd ../../../reconstruction/
% Map the index of the modified Kcat values to the new model (after rxns removals)
modifications = mapModifiedRxns(modifications,model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modified = mapModifiedRxns(modifications,model)
modified = [];
for i=1:length(modifications{1})
    rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
    str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
    modified = [modified; str];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the top growth limiting enzymes that were detected by the
% modifyKcats.m script in a preliminary run.
function [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications)
newValue = [];
reaction = string(reaction);
% chorismate synthase (Q9KXQ4/EC4.2.3.5) - model specified substrate 
% 5-O-(1-Carboxyvinyl)-3-phosphoshikimate is synonymous to the BRENDA
% specified substrate 5-enolpyruvylshikimate 3-phosphate. Changed manually
% to match kcat of wild-type N. crassa enyzme (doi: 10.1111/
% j.1742-4658.2008.06305.x)
if strcmpi('prot_Q9KXQ4',enzName)
    if contains(reaction,'chorismate synthase (')
        newValue         = -(0.87*3600)^-1; % BRENDA: WT N. crassa
        modifications{1} = [modifications{1}; 'Q9KXQ4'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% phosphoribosylformylglycinamidine synthase (Q9RKK5,6,7/EC6.3.5.3) - kcat
% automatically suggested was for ammonium as substrate, not glutamine.
% Instead, use specific activity as provided in BRENDA, which was assayed
% with glutamine according to original paper (doi:10.1021/bi00432a017)
if contains(reaction,'phosphoribosylformylglycinamidine synthase (')
    if strcmpi('prot_Q9RKK5',enzName)
        newValue         = -(5.05*3600)^-1; % BRENDA: WT E. coli
        modifications{1} = [modifications{1}; 'Q9RKK5'];
        modifications{2} = [modifications{2}; reaction];
    end
    if strcmpi('prot_Q9RKK6',enzName)
        newValue         = -(5.05*3600)^-1; % BRENDA: WT E. coli
        modifications{1} = [modifications{1}; 'Q9RKK6'];
        modifications{2} = [modifications{2}; reaction];
    end
    if strcmpi('prot_Q9RKK7',enzName)
        newValue         = -(5.05*3600)^-1; % BRENDA: WT E. coli
        modifications{1} = [modifications{1}; 'Q9RKK7'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% methylmalonate-semialdehyde dehydrogenase (malonic semialdehyde)
% (Q9L1J1/EC1.2.1.27) - kcat automatically assigned was from archaea,
% instead use value from Bacillus subtilis (doi:10.1074/jbc.M110.213280).
if strcmpi('prot_Q9L1J1',enzName)
    if contains(reaction,'methylmalonate-semialdehyde dehydrogenase (malonic semialdehyde) (')
        newValue         = -(2.2*3600)^-1; % BRENDA: WT B. subtilis
        modifications{1} = [modifications{1}; 'Q9L1J1'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% phosphoribosyl-ATP pyrophosphatase (Q9EWK0/EC3.6.1.31) - kcat
% automatically assigned was calculated from specific activity in
% Salmonella enterica, but the reported value was measured in cell extract,
% not from purified enzyme. Instead, use specific activity from
% S. cerevisiae (PMID:379004).
if strcmpi('prot_Q9EWK0',enzName)
    if contains(reaction,'phosphoribosyl-ATP pyrophosphatase (')
        newValue         = -(526*3600)^-1; % BRENDA: WT E. coli
        modifications{1} = [modifications{1}; 'Q9EWK0'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% glyceraldehyde-3-phosphate dehydrogenase (Q9Z518/EC1.2.1.12) - assigned
% kcat from Corynebacterium glutamicum was highly growth limiting.
% Instead use specific activity measured of pentalenolactone sensitive
% gapdh in Streptomyces arenae (PMID:6822480)
if strcmpi('prot_Q9Z518',enzName)
    if contains(reaction,'glyceraldehyde-3-phosphate dehydrogenase (')
        newValue         = -(112*60*MW_set)^-1;
        modifications{1} = [modifications{1}; 'Q9Z518'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% phosphoserine phosphatase (L-serine) (Q9S281/EC3.1.3.3) - assigned
% kcat from Mycobacterium tuberculosis was highly growth limiting.
% Match with pseudonym of substrate L-phosphoserine, use reported kcat
% from Porphyromonas gingivalis (1508 min-1) (PMID:16832066)
if strcmpi('prot_Q9S281',enzName)
    if contains(reaction,'phosphoserine phosphatase (L-serine) (')
        newValue         = -(1508*60)^-1;
        modifications{1} = [modifications{1}; 'Q9S281'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% thymidylate synthase (Flavin-dependent) (O86840/EC2.1.1.148) -
% assigned kcat from Helicobacter pylori was growth limiting. Recent
% new kcat value from E. coli (not yet in BRENDA) (PMID:29715278)
if strcmpi('prot_O86840',enzName)
    if contains(reaction,'thymidylate synthase (Flavin-dependent) (')
        newValue         = -(8.7*3600)^-1;
        modifications{1} = [modifications{1}; 'O86840'];
        modifications{2} = [modifications{2}; reaction];
    end
end
% ATP synthase (Q9K4D5,Q9K4D4,P0A300/EC3.6.3.14) - assigned activity was 
% growth limiting. Use kcat from Bacillus sp. (PMID:24876384)
if (strcmpi('prot_Q9K4D5',enzName) || strcmpi('prot_P0A300',enzName) || ...
    strcmpi('prot_Q9K4D4',enzName)) && ...
    contains(reaction,'ATP synthase (four protons for one ATP) (')
    newValue      = -(390*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9K4D5'];
    modifications{2} = [modifications{2}; reaction];
end
% cytochrome oxidase bd (menaquinol-9: 2 protons) (Q9K451/EC1.9.3.1) -
% assigned activity was growth limiting. Take highest reported kcat for
% this EC number
if strcmpi('prot_Q9K451',enzName) && ...
    contains(reaction,'cytochrome')
    newValue      = -(2000*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9K451'];
    modifications{2} = [modifications{2}; reaction];
end
% 3-oxoacyl-[acyl-carrier-protein] synthase (Q9Z4Y3/EC2.3.1.179) -
% assigned activity from E. coli was growth limiting, 37% of total
% protein in pool-model. Use trimmed mean (10%) of all BRENDA kcat values
% for EC2.3.1.-, which is 33.54.
% Alternatively use S.A., 7.5 umol/min/mg protein and 
% 76 kDA, as reported in PMID: 237914.
% Alternative, not in BRENDA, but activity was also measured in 
% Streptomyces coelicolor (PMID:22136753). Different sets of substrates are
% assayed in that paper, take the highest reported activity of 20.33 min-1.
%%if (strcmpi('prot_Q9Z4Y3',enzName) | strcmpi('prot_Q9K3H4',enzName) | ...
   %% strcmpi('prot_Q9RDP7',enzName) | strcmpi('prot_Q9RK62',enzName)) && ...
if    (contains(reaction,'3-oxoacyl-[acyl-carrier-protein] synthase (') | ...
    contains(reaction,'beta-ketoacyl-ACP synthase (') | ...
    contains(reaction,'3-Oxo-glutaryl-[ACP] methyl ester synthase') | ...
    contains(reaction,'3-Oxo-pimeloyl-[ACP] methyl ester synthase'))    
    newValue      = -(2000*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9Z4Y3'];
    modifications{2} = [modifications{2}; reaction];
end
% demethylmenaquinone methyltransferase (Q9XAP8/EC2.1.1.163) - assigned
% activity was growth limiting. No activities measured. Take the highest
% kcat reported for Streptomyces spp with EC2.1.1.-, this is 8.15 for
% phenylpyruvate C3-methyltransferase (EC2.1.1.281) in Streptomyces
% hygroscopicus (PMID:19731276)
if strcmpi('prot_Q9XAP8',enzName) && ...
    (contains(reaction,'S-adenosylmethione:2-demthylmenaquinole methyltransferase (menaquinone 9)') | ...
    contains(reaction,'demethylmenaquinone methyltransferase'))
    newValue      = -(8.15*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9XAP8'];
    modifications{2} = [modifications{2}; reaction];
end
% GTP cyclohydrolase 1 (Q9X8I3/EC3.5.4.16) - one of the topUsedEnzymes.
% Use S.A. measured from Streptomyces tubercidicus (PMID:9868539), 9.1
% umol.
if strcmpi('prot_Q9X8I3',enzName) && contains(reaction,'GTP cyclohydrolase I')
    newValue      = -(9.1*60*MW_set)^-1;
    modifications{1} = [modifications{1}; 'Q9X8I3'];
    modifications{2} = [modifications{2}; reaction];
end
% geranylgeranyltranstransferase (Q9F2X8/EC2.5.1.31) - no reliable kcat or
% SA (only kcat known is from a plant). Rather, use same Kcat that is
% estimated for other reactions catalyzed by this enzyme: 18.9998
if strcmpi('prot_Q9F2X8',enzName) && contains(reaction,'geranylgeranyltranstransferase')
    newValue      = -(18.9998*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9F2X8'];
    modifications{2} = [modifications{2}; reaction];
end
% Undecaprenyl diphosphate synthase (Q9L2H4/EC2.5.1.31) - EC number on UniProt has
% wild-card and was therefore not queried for kcat value. Take highest
% bacterial kcat value from BRENDA: 2.5 sec-1 for E. coli.
if strcmpi('prot_Q9L2H4',enzName) && contains(reaction,'Undecaprenyl diphosphate synthase')
    newValue      = -(2.5*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9L2H4'];
    modifications{2} = [modifications{2}; reaction];
end
% cardiolipin synthase (Q9KZP3/EC2.7.8.41) - one of the topUsedEnzymes.
% Specific activity was also measured in E. coli (not included in BRENDA)
% as 29 nmol/min/mg protein (PMID:8262233). Assume 54.8 kDa molecular
% weight, k-cat is 0.026487 sec-1
if strcmpi('prot_Q9KZP3',enzName) && ...
        contains(reaction,'cardiolipin synthase II')
    newValue      = -(0.026487*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9KZP3'];
    modifications{2} = [modifications{2}; reaction];
end
% 4-hydroxy-2,2'-bipyrrole-5-methanol methyltransferase (O54158/EC2.1.1.-) 
% - dominates during RED production, kcat is not from a measured value.
% Take the highest kcat reported for Streptomyces spp with EC2.1.1.-, this
% is 8.15 for phenylpyruvate C3-methyltransferase (EC2.1.1.281) in
% Streptomyces hygroscopicus (PMID:19731276)
if strcmpi('prot_O54158',enzName) && ...
        contains(reaction,'4-hydroxy-2,2''-bipyrrole-5-methanol methyltransferase')
    newValue      = -(8.15*3600)^-1;
    modifications{1} = [modifications{1}; 'O54158'];
    modifications{2} = [modifications{2}; reaction];
end
% L-prolyladenylate synthase (O54154/EC6.2.1.53) - assigned kcat values too
% low to support RED production when proteomics data integrated. No kcat
% values are provided in BRENDA, take trimmed mean (10%) of EC6.2.1.-
% as reported on BRENDA: 23.01 sec-1, instead of initially assigned 3.2.
if strcmpi('prot_O54154',enzName) && ...
        contains(reaction,'L-prolyladenylate synthase')
    newValue      = -(23.01*3600)^-1;
    modifications{1} = [modifications{1}; 'O54154'];
    modifications{2} = [modifications{2}; reaction];
end
% cobaltochelatase (Q9RJ19/EC6.6.1.2) - rate limiting. No measured
% activity, take highest report activity of EC6.6.1.-: 0.03
if strcmpi('prot_Q9RJ19',enzName) && ...
        contains(reaction,'cobaltochelatase')
    newValue      = -(0.03*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9RJ19'];
    modifications{2} = [modifications{2}; reaction];
end
% Phosphoadenosine phosphosulfate reductase (Q9ADG3/EC1.8.4.8) - rate limiting. No measured
% activity, take highest report activity of EC6.6.1.-: 0.03
if strcmpi('prot_Q9ADG3',enzName) && ...
        contains(reaction,'phosphoadenylyl-sulfate reductase')
    newValue      = -(3.5*3600)^-1;
    modifications{1} = [modifications{1}; 'Q9ADG3'];
    modifications{2} = [modifications{2}; reaction];
end
% dUTP diphosphatase (No1) (O54134/EC3.6.1.23) - rate limiting, measured at 20 C. 
% Take the highest prokaryotic activity measured at 30 C: 38.3 sec-1 as kcat
% for E. coli, as reported in PMID:346589
if strcmpi('prot_O54134',enzName) && ...
        contains(reaction,'dUTP diphosphatase')
    newValue      = -(38.3*3600)^-1;
    modifications{1} = [modifications{1}; 'O54134'];
    modifications{2} = [modifications{2}; reaction];
end
end
