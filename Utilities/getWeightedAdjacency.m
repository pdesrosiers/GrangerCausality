function [ weightMat ] = getWeightedAdjacency( G )
%%
% INPUT:    G           graph or digraph object
% OUTPUT:   weightMat   weighted adjacency matrix is sparse form

%% Find edges' endpoints and weight
    from = G.Edges.EndNodes(:,1);
    to = G.Edges.EndNodes(:,2);
    weight = G.Edges.Weight;

%% Build the sparse version of the weight matrix
    if isa(G, 'graph')  % if graph, edges muts be symmetrized
        nondiag = from~=to;
        newfrom = from(nondiag);
        newto = to(nondiag);
        newweight = weight(nondiag);    
        from = [ from ; newto];
        to = [ to; newfrom];
        weight = [ weight ; newweight];
    end
    weightMat = sparse(from,to, weight);

end
