function TestCRAveraging()


    mesh = Mesh.loadFromGeometry("unitsquare");
    mesh.refineUniform(1);
    % fesS1 = FeSpace(mesh, LowestOrderH1Fe());
    fesCR = FeSpace(mesh, LowestOrderCRFe());
    % P = LoMeshProlongation(fesCR);
    P = LoMeshAveraging(fesCR);

    figure(3); plot(mesh, 'labelEdges', true);


    uCR = FeFunction(fesCR);
    uCR.setData(0);
    freeDofsCR = getFreeDofs(fesCR);
    nFreeDofsCR = length(freeDofsCR);
    tmp = zeros(nFreeDofsCR, 1);
    idx = randi(nFreeDofsCR, 1);
    % freeDofsCR(idx)
    tmp(idx) = 1;
    uCR.setFreeData(tmp);

    figure(1);
    plot(uCR);

    idxElem = mesh.edge2elements(:,freeDofsCR(idx));
    mesh.refineLocally(idxElem(idxElem~=0));

    uCRFine = FeFunction(fesCR);
    uCRFine.setData(0);

    tmp = prolongate(P, uCR);
    uCRFine.setData(tmp);


for j = 1:5

    mesh.refineUniform();

    tmp = prolongate(P, uCRFine);
    uCRFine.setData(tmp);

end

figure(2);
plot(uCRFine);

end