function dofs = assemblePatchDofs(obj,p)

%% preallocate arrays
mesh = obj.mesh;
dofs = getDofs(obj);
freedofs = getFreeDofs(obj);
patchDofs = {};



for ver = 1:mesh.nCoordinates
    elemPatch = ceil(find(mesh.elements(:)==ver)./3); 
    patchDofs_temp = unique(dofs.element2Dofs(:, elemPatch)); %all dofs: need to exlude boundary

    patchDofs_temp = intersect(patchDofs_temp, freedofs); %only keeping freenodes

    verPatch = unique(mesh.elements(:,elemPatch));
    otherverPatch = setdiff(verPatch, ver);
    patchDofs_temp = setdiff(patchDofs_temp, otherverPatch); %only one vertex contribution -> from ver

    %only edge contributions from edges containing vertex ver
    id = [];
    for neigh = otherverPatch'
        id = [id;  ceil(find(dofs.edge2Dofs(:)==neigh)./(p+1))];
    end
    id = setdiff(unique(id),ceil(find(dofs.edge2Dofs(:)==ver)./(p+1)));
    patchDofs{ver} = setdiff(patchDofs_temp,unique(dofs.edge2Dofs(:,id)));

end


dofs.patch2Dofs = patchDofs; %DOFs on patches

end
