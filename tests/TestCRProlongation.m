function TestCRProlongation()


    mesh = Mesh.loadFromGeometry("unitsquare");
    fesS1 = FeSpace(mesh, LowestOrderH1Fe());
    fesCR = FeSpace(mesh, LowestOrderCRFe());
    P = LoMeshProlongation(fesCR);

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

    uS1 = FeFunction(fesS1);
    uS1.setData(0);


    tmp = prolongate(P, uCR);
    uS1.setData(tmp);

    figure(2);
    plot(uS1);


end