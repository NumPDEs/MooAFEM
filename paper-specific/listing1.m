% ******************************************************************************
% Simple AFEM algorithm for the Poisson equation (Listing 1 in the
% accompanying paper).
% ******************************************************************************

mesh = Mesh.loadFromGeometry('unitsquare');
fes = FeSpace(mesh, LowestOrderH1Fe);
u = FeFunction(fes);
blf = BilinearForm(fes);
blf.a = Constant(mesh, 1);
lf = LinearForm(fes);
lf.f = Constant(mesh, 1);

while mesh.nElements < 1e6
	A = assemble(blf);
	F = assemble(lf);
	freeDofs = getFreeDofs(fes);
	u.setFreeData(A(freeDofs,freeDofs) \ F(freeDofs));
	
	hT = sqrt(getAffineTransformation(mesh).area);
	qrEdge = QuadratureRule.ofOrder(1, '1D');
	qrTri = QuadratureRule.ofOrder(1, '2D');
	volumeRes = integrateElement(CompositeFunction(@(x) x.^2, lf.f), qrTri);
	edgeRes = integrateNormalJump(Gradient(u), qrEdge, @(j) j.^2, {}, ':');
	edgeRes(mesh.boundaries{:}) = 0;
	eta2 = hT.^2.*volumeRes + hT.*sum(edgeRes(mesh.element2edges), Dim.Vector);
	
	marked = markDoerflerSorting(eta2, 0.5);
	mesh.refineLocally(marked, 'NVB');
end