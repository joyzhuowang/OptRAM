function [simplemodel] =  pre_process(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The pre-process of optimization algorithms
% INPUTS
%model
%
%OUTPUTS
%A simple model, remove reactions with lb=ub=0
%

[lb,ub] = fluxVariability(model,0);
t=lb+ub+lb.*ub;%Indicator, t=0 when lb=ub=0
%remove elements for reactions
rmMAT=model.S(:,t~=0);
simplemodel.rxns=model.rxns(t~=0);
if isfield(model,'rev')
    simplemodel.rev=model.rev(t~=0);
end
simplemodel.lb=model.lb(t~=0);
simplemodel.ub=model.ub(t~=0);
simplemodel.c=model.c(t~=0);
simplemodel.grRules=model.grRules(t~=0);
simplemodel.rxnNames=model.rxnNames(t~=0);
%remove elements for metabolites (corresponding to deleted reactions)
x=sum(abs(rmMAT),2);
simplemodel.S=rmMAT(x~=0,:);
simplemodel.mets=model.mets(x~=0);
simplemodel.metNames=model.metNames(x~=0);
simplemodel.b=model.b(x~=0);
%remove corresponding genes
rgMAT=model.rxnGeneMat(t~=0,:);
y=sum(abs(rgMAT));
simplemodel.rxnGeneMat=rgMAT(:,y~=0);
simplemodel.genes=model.genes(y~=0);

end