function [simpleindmax,fit] =  post_process(indmax,new_model,v,regnet,BPCYids)
%Post-processing, to remove sites that do not affect score of objective function
%
% INPUTS
%indmax--original final solution
%model--
%v--
%regnet--
%BPCYids--
%
% OUTPUTS
% simpleindmax--Simplified final solution
%

if numel(indmax.gene)==0
    simpleindmax=0;
    fit=0;
else
[objmax,~,~,target,growth,~,~] = OptRAM_obj(indmax,new_model,v,regnet,BPCYids);
x1=round(objmax*10000)/10000;
n=numel(indmax.gene);
fit=zeros(n+1,3);
delid=[];
fit(1,1)=x1;
fit(1,2)=target;
fit(1,3)=growth;
for i=1:n
    tempindmax=indmax;
    tempindmax.code(i)=6;
    [tempfit,~,~,target,growth,~,~] = OptRAM_obj(tempindmax,new_model,v,regnet,BPCYids);
    fit(i+1,1)=tempfit;
    fit(i+1,2)=target;
    fit(i+1,3)=growth;
    x2=round(tempfit*10000)/10000;
    if x2>=x1*0.99
        delid=[delid i];
    end
end 
simpleindmax=indmax;
simpleindmax.gene(delid)=[];
simpleindmax.code(delid)=[];
fit(delid,:)=[];
end
end