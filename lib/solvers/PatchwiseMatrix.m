classdef PatchwiseMatrix < handle
    properties(Access=protected)
        fes
        patch
        patchElements
        ptr
        preA
        matrix
    end

    methods(Access=public)
        function obj = PatchwiseMatrix(fes, preA);
            obj.fes = fes;
            [obj.patch, obj.patchElements, obj.ptr] = assemblePatchDofs(obj.fes);
            obj.preA = preA;

            obj.assemble();
        end
    end

    methods(Access=protected)
        function assemble(obj)
            
            D = getDofs(obj.fes).element2Dofs;
            ElPatches = obj.patchElements;
            ptr = obj.ptr;


            N = getDofs(obj.fes).nDofs;
            n = length(D(1,:)); %number of elements
            Nc = length(ptr)-1; %number of vertices
            ndofEl = length(D(:,1)); %number of dofs per element

            %patchwiseCholmodif = cell(Nc,1);

            patch_sizes = cellfun(@numel, obj.patch);
            tot_patchsizes = sum(patch_sizes.^2);
            I = zeros(tot_patchsizes,1);
            J = zeros(tot_patchsizes,1);

            idx = 0;
            global2local = zeros(getDofs(obj.fes).nDofs, 1);
            for i = 1:Nc
                
                idx = idx(end) + (1:(patch_sizes(i)^2));

                idxD = obj.patch{i};
                global2local(idxD) = 1:patch_sizes(i);
                idxE = ElPatches(ptr(i):ptr(i+1)-1);
                localIdx = reshape(global2local(D(:,idxE)), ndofEl, []);
                allDofs = D(:,idxE); %note: NOT freedof numbering
                freeDofspatch = obj.patch{i}; %which global freeDofs on patch i

                I(idx) = repmat(localIdx, 1, patch_sizes(i));
                J(idx) = repelem(localIdx, 1, patch_sizes(i));

                patchwiseA
            end

            localA = sparse(I, J, patchwiseA);
            patchwiseCholmodif = chol(localA);
        end


    end
end

