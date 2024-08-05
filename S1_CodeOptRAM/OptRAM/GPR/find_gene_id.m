function [ position ] = find_gene_id( model , gene_name )

position = -1;
for i = 1 : 1 : length(model.genes)
    if strcmp(gene_name , model.genes(i))
        position = i;
        return;
    end
end
end