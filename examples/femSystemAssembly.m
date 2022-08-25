% ******************************************************************************
% Example that showcases how to assemble FEM linear systems.
% ******************************************************************************

%% initialisation
% First, one needs a mesh and a finite element space.
mesh = Mesh.loadFromGeometry('unitsquare');
mesh.refineUniform(4);
fes = FeSpace(mesh, LowestOrderH1Fe());

% One can define a bilinear form 'A' and a linear form 'F', such that
%   A(u,v) = \int a*Du*Dv + b*Du*v + c*u*v dx,
%   F(v) = \int f*v + fvec*Dv dx.
A = BilinearForm();
F = LinearForm();

%% coefficients
% The coefficients a, b, c and f, fvec can be arbitrary Evaluables.
A.a = Constant(mesh, [1;0;0;10]);
F.f = Constant(mesh, 1);

%% assembly and solution
% Once the coefficients are fixed, the bilinear form and linear form can be
% assembled to yield a (sparse) matrix and a vector, respectively.
Amat = assemble(A, fes);
Fvec = assemble(F, fes);

% This allows for easy solution of the corresponding linear system.
freeDofs = getFreeDofs(fes);
u = FeFunction(fes);
u.setFreeData( Amat(freeDofs, freeDofs) \ Fvec(freeDofs) );
plot(u)

%% change of coefficients
% Coefficients can be set/unset at any time. They are evaluated only at assembly.
A.a = Constant(mesh, 1);
A.b = Gradient(u);
A.c = CompositeFunction(@(x) x.^2, u);
F.f = [];
F.fvec = MeshFunction(mesh, @(x) -(0.5-x).^2);

Amat = assemble(A, fes);
Fvec = assemble(F, fes);

w = FeFunction(fes);
w.setFreeData( Amat(freeDofs, freeDofs) \ Fvec(freeDofs) );
plot(w)

%% other FeSpaces
% Of course, also other FeSpaces can be used to assemble the forms.
fes2 = FeSpace(mesh, LowestOrderL2Fe());
M = BilinearForm();
G = LinearForm();

% This realizes mass matrix and integration for piecewise constant functions.
% NOTE: All assembly operations are implemented in terms of local shape
% functions. There is no check whether certain coefficients (e.g. non-vanishing
% diffusion for discontinuous, piecewise constant elements) are sensible.
M.c = Constant(mesh, 1);
G.f = Constant(mesh, 1);

Mmat = assemble(M, fes2);
Gvec = assemble(G, fes2);

v = FeFunction(fes2);
v.setData(ones(mesh.nElements,1));

fprintf('F(1) = %.3f = %.3f = M(1,1)\n', v.data*Gvec, v.data*Mmat*v.data')
