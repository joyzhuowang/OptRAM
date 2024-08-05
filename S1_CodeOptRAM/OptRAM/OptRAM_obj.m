function [obj,rxn_p,mgene_p,target,growth,FBAsolution,model2] = OptRAM_obj(individual,model,v,regnet,BPCYids)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The objective function of optimization algorithms
% INPUTS
%individual--Two-dimensional array of mutant, including genes and codes
%model--Metabolic model
%v--Reference flux value
%regnet--Regulatory network
%BPCYids--Reaction IDs in BPCY
%
% OUTPUTS
% obj--Objective function to be maximized 
% target--Flux value of target product reaction
% growth--Flux value of biomass
% model2--New metabolic model of the mutant
%
mgene_p=ones(numel(model.genes),1);
rules=model.grRules;
expression=model.postfix_expression;
long=model.postfix_long;
n_mgene=numel(model.genes);
for i=1:length(individual.gene)
    if individual.gene(i) <= n_mgene
        mgene_p(individual.gene(i))=precode(individual.code(i));
    end
end
if ~isempty(individual.gene(individual.gene>numel(regnet.mgene)))
    mgene_p=gene_cal(individual,regnet,mgene_p,n_mgene);
end
n=length(model.rxns);
rxn_p=ones(n,1);
for i=1:n
    if ~isempty(rules{i})
        temp = calc_postfix(expression,long,i,mgene_p);
        rxn_p(i) = temp(1);
    end
end
%Change constraints according to reaction expression changes
for i=1:n
    p=rxn_p(i);
    if p>1
        if v(i)>0
            model.lb(i)=min(p*v(i),model.ub(i));
        elseif v(i)<0
            model.ub(i)=max(p*v(i),model.lb(i));
        end
    elseif p<1
        if v(i)>0
            model.ub(i)=p*v(i);
            model.lb(i)=0;
        elseif v(i)<0
            model.lb(i)=p*v(i);
            model.ub(i)=0;
        else
            model.lb(i)=0;
            model.ub(i)=0;
        end
    end
end
%model = changeObjective(model,model.rxns(BPCYids(1)));
FBAsolution = optimizeCbModel(model);
model2 = changeObjective(model,model.rxns(BPCYids(2)));
if numel(FBAsolution.x)==0
    obj=0;
    growth=0;
    target=0;
else
    model.lb(BPCYids(1)) = FBAsolution.f*0.99;
    model2 = changeObjective(model,model.rxns(BPCYids(2)));
    FBAsolution1 = optimizeCbModel(model2,'min');
    FBAsolution2 = optimizeCbModel(model2,'max');
    if numel(FBAsolution1.x)==numel(FBAsolution2.x) && numel(FBAsolution1.x)~=0
        FBAsolution.x = (FBAsolution1.x+FBAsolution2.x)./2;
        rang = (FBAsolution2.x(BPCYids(2))-FBAsolution1.x(BPCYids(2)))/2;
    else
        FBAsolution.x = [];
        rang=100;
    end   
    if numel(FBAsolution.x)==0
        obj=0;
        growth=0;
        target=0;
    else
        growth=FBAsolution.x(BPCYids(1));%reaction id of growth
        target=FBAsolution.x(BPCYids(2));%reaction id of target
        substrate=FBAsolution.x(BPCYids(3));%reaction id of substrate
        if rang==0
            rang=rang+1E-3;
        end
        obj=abs(growth*target/substrate)*(1-log(rang/target));
    end
end

end


function p=precode(code)
if code==0
    p=0.001;
elseif code==6
    p=1;
else
    p=2^code;
end
end

function mgene_p=gene_cal(individual1,regnet,mgene_p,n_mgene)
%%Calculate fold changes in gene expression based on IDREAM
    TFs=individual1.gene(individual1.gene>n_mgene);
    TFs_code=individual1.code(individual1.gene>n_mgene);
    TFs=TFs-n_mgene;
    for i=1:numel(TFs)
        p=precode(TFs_code(i));
        regnet.TFexp(TFs(i))=regnet.TFexp(TFs(i))*p;
    end
    a=regnet.mat(:,TFs);
    gene2cal=find(sum(abs(a),2)>0);
    calmat=regnet.mat(gene2cal,:);
    p_gene2cal=2.^(calmat*log2(regnet.TFexp));
    mgeneid=regnet.mgeneid;
    mgene_p(mgeneid(gene2cal))=p_gene2cal;
end

