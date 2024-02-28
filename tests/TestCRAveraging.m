function TestCRAveraging()


    mesh = Mesh.loadFromGeometry("unitsquare");
    % fesS1 = FeSpace(mesh, LowestOrderH1Fe());
    fesCR = FeSpace(mesh, LowestOrderCRFe());
    P = LoMeshAveraging(fesCR);

    figure(3); plot(mesh, 'labelEdges', true);

    uCR = FeFunction(fesCR);
    uCR.setData(0);
    freeDofsCR = getFreeDofs(fesCR);
    nFreeDofsCR = length(freeDofsCR);
    tmp = zeros(nFreeDofsCR, 1);
    tmp(freeDofsCR(1)) = 1;
    uCR.setFreeData(tmp);


    figure(1);
    plot(uCR);

    mesh.refineUniform(1);

    uCRFine = FeFunction(fesCR);
    uCRFine.setData(0);


    tmp = prolongate(P, uCR);
    uCRFine.setData(tmp);

    figure(2);
    plot(uCRFine);


end