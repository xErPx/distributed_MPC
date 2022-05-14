function [ partition ] = get_explicit_partition( mpsol )

% Extract partition
mpsol_Partition = mpsol{1}.Pn;

% Number of regions in partition
Nr = length(mpsol{1}.Fi);

%% Construct regions using function Polyhedron
for r = 1 : Nr 
    % Convert MPT2 Polytope into MPT3 Polyhedron
    [H,K] = double(mpsol{1}.Pn(r)); 
    P(r) = Polyhedron('A',H,'b',K);
end % for r

%% Assign outputs
mpsol{1}.Partition = P;
partition = mpsol;

end % function