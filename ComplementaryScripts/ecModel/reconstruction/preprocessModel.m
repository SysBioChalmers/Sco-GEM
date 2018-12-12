%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,name,version] = preprocessModel(model,name,version)
% Performs some preliminary modifications to the metabolic model & 
% retrieves the model's name & version (either by parsing model.id or by
% asking the user to input it), if they were not already defined.
%
% model     A genome-scale model in RAVEN format
% name      The name of the model (alternatively, an empty string)
% version   The version of the model (alternatively, an empty string)
% 
% model     The processed model
% name      The resulting name of the model (if not specified before)
% version   The resulting version of the model (if not specified before)
%
% Benjamin J. Sanchez. Last edited: 2018-09-01
% Eduard Kerkhoven. Last edited: 2018-10-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,name,version] = preprocessModel(model,name,version)

%Remove gene rules from pseudoreactions (if any):
for i = 1:length(model.rxns)
    if endsWith(model.rxnNames{i},' pseudoreaction')
        model.grRules{i}      = '';
        model.rxnGeneMat(i,:) = zeros(1,length(model.genes));
    end
end

%Swap direction of only reactions that are defined to only carry negative flux
to_swap=model.lb < 0 & model.ub == 0;
model.S(:,to_swap)=-model.S(:,to_swap);
model.ub(to_swap)=-model.lb(to_swap);
model.lb(to_swap)=0;

%Delete blocked rxns (LB = UB = 0):
to_remove = logical((model.lb == 0).*(model.ub == 0));
model     = removeReactions(model,model.rxns(to_remove),true,true,true);

%Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
        model.rev(i) = true;
    end
end

%Correct complex grRules

model = changeGeneAssoc(model, 'OXDHCOAT', '(SCO5144 and SCO6701) or (SCO5144 and SCO6967) or (SCO5144 and SCO3079) or (SCO5144 and SCO6731)');
model = changeGeneAssoc(model, 'IBTMr', '(SCO5415 and SCO4800) or (SCO5415 and SCO6833)');
model = changeGeneAssoc(model, 'AKGDH', '(SCO5281 and SCO2181 and SCO0884) or (SCO5281 and SCO2181 and SCO2180) or (SCO5281 and SCO2181 and SCO4919) or (SCO5281 and SCO7123 and SCO0884) or (SCO5281 and SCO7123 and SCO2180) or (SCO5281 and SCO7123 and SCO4919)');
model = changeGeneAssoc(model, 'AKGDH2', '(SCO4594 and SCO4595 and SCO0681) or (SCO6269 and SCO6270 and SCO0681)');
model = changeGeneAssoc(model, 'PAPSR', '(SCO6100 and SCO0885) or (SCO6100 and SCO3889) or (SCO6100 and SCO5419) or (SCO6100 and SCO5438)');
model = changeGeneAssoc(model, 'GLYCL', '(SCO5471 and SCO1378 and SCO5472 and SCO0884) or (SCO5471 and SCO1378 and SCO5472 and SCO2180) or (SCO5471 and SCO1378 and SCO5472 and SCO4919)');
model = changeGeneAssoc(model, '2MBCOATA', '(SCO1271 and SCO2389) or (SCO1271 and SCO0549) or (SCO1271 and SCO1267) or (SCO1271 and SCO1272) or (SCO2388 and SCO2389) or (SCO2388 and SCO0549) or (SCO2388 and SCO1267) or (SCO2388 and SCO1272) or (SCO6564 and SCO2389) or (SCO6564 and SCO0549) or (SCO6564 and SCO1267) or (SCO6564 and SCO1272)');
model = changeGeneAssoc(model, {'ACCOAC','ACCOAC_1'}, '(SCO2445 and SCO2777) or (SCO2445 and SCO4921) or (SCO2445 and SCO6271) or (SCO5535 and SCO5536 and SCO2777) or (SCO5535 and SCO5536 and SCO4921) or (SCO5535 and SCO5536 and SCO6271)');
model = changeGeneAssoc(model, {'ACOATA', 'BCOATA', 'IBCOATA', 'IVCOATA', 'PCOATA'}, '(SCO1271 and SCO2389) or (SCO1271 and SCO0549) or (SCO1271 and SCO1267) or (SCO1271 and SCO1272) or (SCO2388 and SCO2389) or (SCO2388 and SCO0549) or (SCO2388 and SCO1267) or (SCO2388 and SCO1272) or (SCO6564 and SCO2389) or (SCO6564 and SCO0549) or (SCO6564 and SCO1267) or (SCO6564 and SCO1272)');
model = changeGeneAssoc(model, 'MCOATA', '(SCO2387 and SCO2389) or (SCO2387 and SCO0549) or (SCO2387 and SCO1267) or (SCO2387 and SCO1272)');
model = changeGeneAssoc(model, 'METSOXR1', '(SCO4956 and SCO0885) or (SCO4956 and SCO3889) or (SCO4956 and SCO5419) or (SCO4956 and SCO5438)');
model = changeGeneAssoc(model, 'METSOXR2', '(SCO6061 and SCO0885) or (SCO6061 and SCO3889) or (SCO6061 and SCO5419) or (SCO6061 and SCO5438)');
model = changeGeneAssoc(model, {'MPTG','MPTG2'}, '(SCO3847 and SCO2709) or (SCO3847 and SCO3894) or (SCO5301 and SCO2709) or (SCO5301 and SCO3894)');
model = changeGeneAssoc(model, {'RNDR1', 'RNDR2', 'RNDR3', 'RNDR4'},'(SCO5225 and SCO5226 and SCO0885) or (SCO5225 and SCO5226 and SCO3889) or (SCO5225 and SCO5226 and SCO5419) or (SCO5225 and SCO5226 and SCO5438) or (SCO5805 and SCO0885) or (SCO5805 and SCO3889) or (SCO5805 and SCO5419) or (SCO5805 and SCO5438)');
model = changeGeneAssoc(model, 'CYO2a', '(SCO2150 and SCO2149 and SCO7236) or (SCO2150 and SCO2149 and SCO2148) or (SCO2150 and SCO2149 and SCO7120)');
model = changeGeneAssoc(model, 'CYO2b', '(SCO1934 and SCO2156 and SCO2151 and SCO1930 and SCO7234) or (SCO1934 and SCO2156 and SCO2151 and SCO1930 and SCO2155)');
model = changeGeneAssoc(model, 'SUCD3', '(SCO4856 and SCO4855 and SCO4858 and SCO4857) or (SCO4856 and SCO5106 and SCO4858 and SCO4857) or (SCO5107 and SCO4855 and SCO4858 and SCO4857) or (SCO5107 and SCO5106 and SCO4858 and SCO4857) or (SCO7109 and SCO4855 and SCO4858 and SCO4857) or (SCO7109 and SCO5106 and SCO4858 and SCO4857)');
model = changeGeneAssoc(model, 'TRDR', '(SCO3890 and SCO0885) or (SCO3890 and SCO3889) or (SCO3890 and SCO5419) or (SCO3890 and SCO5438) or (SCO6834 and SCO0885) or (SCO6834 and SCO3889) or (SCO6834 and SCO5419) or (SCO6834 and SCO5438) or (SCO7298 and SCO0885) or (SCO7298 and SCO3889) or (SCO7298 and SCO5419) or (SCO7298 and SCO5438)');
model = changeGeneAssoc(model, 'PPCOAC', '(SCO2776 and SCO2777) or (SCO4380 and SCO4381) or (SCO4921 and SCO4925) or (SCO4921 and SCO4926) or (SCO6271 and SCO4925) or (SCO6271 and SCO4926)');
model = changeGeneAssoc(model, '2OXOADOX', '(SCO5281 and SCO2181 and SCO0884) or (SCO5281 and SCO2181 and SCO2180) or (SCO5281 and SCO2181 and SCO4919) or (SCO5281 and SCO7123 and SCO0884) or (SCO5281 and SCO7123 and SCO2180) or (SCO5281 and SCO7123 and SCO4919)');
model = changeGeneAssoc(model, 'THIORDXi', '(SCO2901 and SCO0885) or (SCO7353 and SCO0885) or (SCO2901 and SCO3889) or (SCO7353 and SCO3889) or (SCO2901 and SCO5419) or (SCO7353 and SCO5419) or (SCO2901 and SCO5438) or (SCO7353 and SCO5438) or SCO4444');
model = changeGeneAssoc(model, {'OIVD1r', 'OIVD2', 'OIVD3'}, '(SCO3816 and SCO3817 and SCO3815 and SCO0884) or (SCO3816 and SCO3817 and SCO3815 and SCO2180) or (SCO3816 and SCO3817 and SCO3815 and SCO4919) or (SCO3816 and SCO3817 and SCO3829 and SCO0884) or (SCO3816 and SCO3817 and SCO3829 and SCO2180) or (SCO3816 and SCO3817 and SCO3829 and SCO4919) or (SCO3830 and SCO3831 and SCO3815 and SCO0884) or (SCO3830 and SCO3831 and SCO3815 and SCO2180) or (SCO3830 and SCO3831 and SCO3815 and SCO4919) or (SCO3830 and SCO3831 and SCO3829 and SCO0884) or (SCO3830 and SCO3831 and SCO3829 and SCO2180) or (SCO3830 and SCO3831 and SCO3829 and SCO4919)');
model = changeGeneAssoc(model, {'ACHBS','ACLS'},'(SCO5512 and SCO5513) or (SCO2769 and SCO5513) or (SCO6584 and SCO5513)');
model = changeGeneAssoc(model, 'ATNS_nh4','(SCO3213 and SCO3214) or (SCO3213 and SCO2043)');
model = changeGeneAssoc(model, 'CDAS2','(SCO3246 and SCO3249) or (SCO1271 and SCO3249) or (SCO2388 and SCO3249) or (SCO6564 and SCO3249) or (SCO6826 and SCO3249) or (SCP1233B and SCO3249)');
model = changeGeneAssoc(model, 'CU1O','SCO3439 and SCO3440 and SCO6712');
model = changeGeneAssoc(model, 'CYO1ab','(SCO2150 and SCO2149 and SCO7236) or (SCO2150 and SCO2149 and SCO2148) or (SCO2150 and SCO2149 and SCO7120)');
model = changeGeneAssoc(model, {'GLUSx','GLUSy'},'(SCO1977 and SCO2026) or (SCO2025 and SCO2026)');
model = changeGeneAssoc(model, 'KAS15','(SCO1271 and SCO2389) or (SCO1271 and SCO0549) or (SCO1271 and SCO1267) or (SCO1271 and SCO1272) or (SCO2388 and SCO2389) or (SCO2388 and SCO0549) or (SCO2388 and SCO1267) or (SCO2388 and SCO1272) or (SCO6564 and SCO2389) or (SCO6564 and SCO0549) or (SCO6564 and SCO1267) or (SCO6564 and SCO1272) or (SCO3246 and SCO2389) or (SCO3246 and SCO0549) or (SCO3246 and SCO1267) or (SCO3246 and SCO1272) or (SCO6826 and SCO2389) or (SCO6826 and SCO0549) or (SCO6826 and SCO1267) or (SCO6826 and SCO1272) or (SCP1233B and SCO2389) or (SCP1233B and SCO0549) or (SCP1233B and SCO1267) or (SCP1233B and SCO1272)');
model = changeGeneAssoc(model, 'MMM','(SCO4800 and SCO4869) or (SCO4800 and SCO6832) or (SCO4800 and SCO4515) or (SCO6833 and SCO4869) or (SCO6833 and SCO6832) or (SCO6833 and SCO4515)');
model = changeGeneAssoc(model, 'SUCD9','(SCO0922 and SCO0923) or (SCO0922 and SCO4856) or (SCO0922 and SCO5107) or (SCO0922 and SCO7109) or (SCO4855 and SCO0923) or (SCO4855 and SCO4856) or (SCO4855 and SCO5107) or (SCO4855 and SCO7109) or (SCO5106 and SCO0923) or (SCO5106 and SCO4856) or (SCO5106 and SCO5107) or (SCO5106 and SCO7109)');

% Redefine, subunits were missing
model = changeGeneAssoc(model, 'PDH', '(SCO1269 and SCO1270 and SCO2183 and SCO0884 and SCO1268) or (SCO1269 and SCO1270 and SCO2183 and SCO0884 and SCO7123) or (SCO1269 and SCO1270 and SCO2183 and SCO0884 and SCO2181) or (SCO1269 and SCO1270 and SCO2183 and SCO2180 and SCO1268) or (SCO1269 and SCO1270 and SCO2183 and SCO2180 and SCO7123) or (SCO1269 and SCO1270 and SCO2183 and SCO2180 and SCO2181) or (SCO1269 and SCO1270 and SCO2183 and SCO4919 and SCO1268) or (SCO1269 and SCO1270 and SCO2183 and SCO4919 and SCO7123) or (SCO1269 and SCO1270 and SCO2183 and SCO4919 and SCO2181) or (SCO1269 and SCO1270 and SCO2371 and SCO0884 and SCO1268) or (SCO1269 and SCO1270 and SCO2371 and SCO0884 and SCO7123) or (SCO1269 and SCO1270 and SCO2371 and SCO0884 and SCO2181) or (SCO1269 and SCO1270 and SCO2371 and SCO2180 and SCO1268) or (SCO1269 and SCO1270 and SCO2371 and SCO2180 and SCO7123) or (SCO1269 and SCO1270 and SCO2371 and SCO2180 and SCO2181) or (SCO1269 and SCO1270 and SCO2371 and SCO4919 and SCO1268) or (SCO1269 and SCO1270 and SCO2371 and SCO4919 and SCO7123) or (SCO1269 and SCO1270 and SCO2371 and SCO4919 and SCO2181) or (SCO1269 and SCO1270 and SCO7124 and SCO0884 and SCO1268) or (SCO1269 and SCO1270 and SCO7124 and SCO0884 and SCO7123) or (SCO1269 and SCO1270 and SCO7124 and SCO0884 and SCO2181) or (SCO1269 and SCO1270 and SCO7124 and SCO2180 and SCO1268) or (SCO1269 and SCO1270 and SCO7124 and SCO2180 and SCO7123) or (SCO1269 and SCO1270 and SCO7124 and SCO2180 and SCO2181) or (SCO1269 and SCO1270 and SCO7124 and SCO4919 and SCO1268) or (SCO1269 and SCO1270 and SCO7124 and SCO4919 and SCO7123) or (SCO1269 and SCO1270 and SCO7124 and SCO4919 and SCO2181)');

% Not standardized, due to too many combinations
%model = changeGeneAssoc(model, {'NADH17b','NADH8'}, 'SCO4564 and SCO4566 and SCO4568 and (SCO4562 or SCO4599) and (SCO4563 or SCO4600) and (SCO3392 or SCO4565) and (SCO4567 or SCO6560) and (SCO4569 or SCO4602) and (SCO4570 or SCO4603) and (SCO4571 or SCO4604) and (SCO4572 or SCO4605) and (SCO4573 or SCO4606 or SCO6954) and (SCO4574 or SCO4607) and (SCO4575 or SCO4608)');
model = changeGeneAssoc(model, {'NADH17b','NADH8'}, '(SCO4562 and SCO4563 and SCO4564 and SCO4565 and SCO4566 and SCO4567 and SCO4568 and SCO4569 and SCO4570 and SCO4571 and SCO4572 and SCO4573 and SCO4574 and SCO4575) or (SCO4562 and SCO4563 and SCO4564 and SCO3392 and SCO4566 and SCO4567 and SCO4568 and SCO4569 and SCO4570 and SCO4571 and SCO4572 and SCO4573 and SCO4574 and SCO4575)');


if isempty(name) && isempty(version) && isfield(model,'id')
    try
        id = strsplit(model.id,'_v');
        if length(id) == 2
            name    = id{1};
            name    = ['ec' upper(name(1)) name(2:end)];
            version = id{2};
        end
    catch
        disp('Not possible to parse name & version. Input manually')
    end
end
while isempty(name)
    name = input('Please enter the desired ecModel name: ','s');
end
while isempty(version)
    version = input('Please enter the model version: ','s');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%