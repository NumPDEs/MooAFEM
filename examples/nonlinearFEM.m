% ******************************************************************************
% Example of solving a nonlinear FEM system with an iterative algorithm.
%   - 1/10*\Delta u + u^3 = f,
% where f is chosen such that the exact solution is
%   u(x,y) = sin(pi*x)*sin(2pi*y).
% ******************************************************************************

%% initialization
nMax = 1000;
tol = 1e-3;

mesh = Mesh.loadFromGeometry('unitsquare');
mesh.refineUniform(5, 'RGB');
fes = FeSpace(mesh, LowestOrderH1Fe());
u = FeFunction(fes);
dofs = getDofs(fes);
u.setData( zeros(dofs.nDofs, 1) );
residualNorm = zeros(nMax,1);

uexact = MeshFunction(mesh, @(x) prod(sin([pi;2*pi].*x), Dim.Vector));

%% problem data
% One can define multiple (bi-)linear forms in order to separate contribution
% which change and those which don't.
% The following iteration is based on the approximation
% 	(u^{n+1})^3 ~ 3*(u^{n})^2 * u^{n+1} - 2*(u^{n})^3.
blfFixed = BilinearForm();
blfVariable = BilinearForm();
lfFixed = LinearForm();
lfVariable = LinearForm();

blfFixed.a = Constant(mesh, 1e-1);
lfFixed.f = CompositeFunction(@(uex) 1e-1*5*pi^2*uex + uex.^3, uexact);
lfFixed.qrf = QuadratureRule.ofOrder(4);

blfVariable.c = CompositeFunction(@(u) 3*u.^2, u);
blfVariable.qrc = QuadratureRule.ofOrder(4);
lfVariable.f = CompositeFunction(@(u) 2*u.^3, u);
lfVariable.qrf = QuadratureRule.ofOrder(4);

B = assemble(blfFixed, fes);
G = assemble(lfFixed, fes);

%% iterative solution of non-linear equation
n = 1;
while n <= nMax
    % The variable part of the (bi-)linear form needs to be assembled in every
    % iteration, since the underlying data change. Since FeFunctions are
    % handles, the current state of u is taken into account in the assembly.
    A = B + assemble(blfVariable, fes);
    F = G + assemble(lfVariable, fes);  
    
    % If we update the data of u, also coefficients which depend on u are updated.
    freeDofs = getFreeDofs(fes);
    u.setFreeData( A(freeDofs,freeDofs) \ F(freeDofs) );
    
    % Various norms for stopping criteria can be computed with quadrature.
    res = CompositeFunction(@(uex,u) (uex-u).^2, uexact, u);
    qr = QuadratureRule.ofOrder(3);
    elementL2Norm = integrateElement(res, qr);
    residualNorm(n) = sqrt(sum(elementL2Norm));
    printLogMessage('Step %d, L2-norm of residual: %.2e', n, residualNorm(n));
    
    if residualNorm(n) < tol
        break
    end
    
    n = n+1;
end

%% plot
plot(uexact)
plot(u)
zlim([-1,1])
plot(res)
