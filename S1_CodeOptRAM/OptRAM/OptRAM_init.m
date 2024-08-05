function [v,regnet,model] = OptRAM_init(rawmodel,rawregnet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% The initialization of optimization algorithms
% INPUTS
%rawmodel--original metabolic model
%rawregnet--original regulatory network
%
% OUTPUTS
% v--reference flux distribution from pFBA 
% model--processed metabolic model
% regnet--processed regulatory network
%
model =  pre_process(rawmodel);
model = pre_calc(model);
%pFBA,1--only minimize the sum of the flux through gene-associated fluxes
[~,~,modelIrrev,MiniFlux] = pFBA(model,'geneoption',1,'skipclass',1);
v=zeros(numel(model.rxns),1);
cnt=1;
for i = 1:length(modelIrrev.rxns)-1
  if (modelIrrev.match(i) ~= 0)
    % Reversible reaction
    if (strcmpi(modelIrrev.rxns{i}(end-1:end),'_b'))
      v(cnt) = v(cnt)-MiniFlux.x(i);
      cnt = cnt + 1;
    else
      v(cnt) = MiniFlux.x(i);
    end
  else
    % Non-reversible reaction
    v(cnt) = MiniFlux.x(i);
    cnt = cnt + 1;
  end
end
%regnet
mat=rawregnet.mat;
targets=rawregnet.mgene;
regnet.TFnames=rawregnet.TFnames;
[newtargets,newmat] = find_mat_mgene(model,mat,targets);
regnet.mat=newmat;
regnet.mgene=newtargets;
tfremove=find(sum(abs(regnet.mat))==0);
regnet.TFnames(tfremove)=[];
regnet.mat(:,tfremove)=[];
generemove=find(sum(abs(regnet.mat),2)==0);
regnet.mgene(generemove)=[];
regnet.mat(generemove,:)=[];
n=numel(regnet.TFnames);
regnet.TFexp=ones(n,1);
regnet.mgeneid=zeros(numel(regnet.mgene),1);
for i=1:numel(regnet.mgene)
    x=find(strcmpi(model.genes,regnet.mgene(i)));
    regnet.mgeneid(i)=x;
end
end

function [newtargets,newmat] = find_mat_mgene(model,mat,targets)
n=length(model.genes);
newmat=[];
newtargets=model.genes;
cnt=0;
del=[];
for i=1:n
    idx=find(strcmpi(targets,model.genes(i)), 1);
    if ~isempty(idx)
        cnt=cnt+1;
        newmat(cnt,:)=mat(idx,:);
    else
        del=[del,i];
    end
end
newtargets(del)=[];
end