% ******************************************************************************
% Example that showcases how to generate functions and what can be done with them.
% ******************************************************************************

%% initialisation
% First, one needs a mesh to visualize (and even define some) functions.
mesh = Mesh.loadFromGeometry('unitsquare');
mesh.refineUniform(3);

% Every subclass of abstract class Evaluable can be evaluated and therefore
% viewed as a function. Functions can be scalar or have multiple values.
a = Constant(mesh, pi);   % constant pi function
f = MeshFunction(mesh, @(x) x);   % vectorfield with source at 0

%% evaluation and plotting
% Functions can be evaluated at barycentric coordinates
midpoints = eval(f, [1;1;1]/3);

% This is used internally to plot them. Typical patch options can be given.
plot(a, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'LineWidth', 1.5)

% Non-linear functions can also be plotted on a coarse mesh
% SubdivisionFactor roughly prescribes how many internal elements should at
% least be used for plotting (patch options not supported)
coarseMesh = Mesh.loadFromGeometry('Zshape');
g = MeshFunction(coarseMesh, @(x) sum(x.^2, Dim.Vector));
plot(g, 'SubdivisionFactor', 1000)

%% fem functions
% Also functions on finite element spaces can be viewed as Evaluables. FEM
% functions need to be given explicit data, which usually comes from
% interpolation, or solution of linear systems.
fesH1 = FeSpace(mesh, LowestOrderH1Fe());
u = FeFunction(fesH1);

% Set data on free dofs, where other dofs are not changed or initialized to 0.
u.setFreeData(1);
plot(u)

% Set data explicitely on all dofs (which correspond to coordinates here).
vectorField = [0.5;0.5] - mesh.coordinates;
u.setData(vectorProduct(vectorField, vectorField));
plot(u)

% Also, non-continuous functions can be created and plotted.
fesL2 = FeSpace(mesh, LowestOrderL2Fe());
v = FeFunction(fesL2);
v.setData(midpoints(1,:) > 0.5);
plot(v)

%% interpolation
% Functions can be interpolated into finite element spaces (also of higher
% order):
fesH1p5 = FeSpace(coarseMesh, HigherOrderH1Fe(5));
f = MeshFunction(coarseMesh, @(x) prod(sin(pi*x),Dim.Vector));
w = FeFunction(fesH1p5);
w.setData(nodalInterpolation(f, fesH1p5));
plot(w, 'SubdivisionFactor', 10000)

%% derivation
% For FEM functions, the (piecewise) gradient can be defined naturally. It can
% be accessed via the Gradient class.
Du = Gradient(u);
plot(Du)

%% function composition
% The CompositeFunction class allows for convenient composition of functions. It
% does not matter if the composision is linear, or the functions of differnt
% type, which makes it a powerful tool for coefficients of PDEs.
% The following function realizes u * v * |Du|
f = CompositeFunction(@(a,b,c) a .* b .* sqrt(vectorProduct(c, c)), u, v, Du);
plot(f);

%% integration
% Every Evaluable can be integrated (element-wise) over the mesh, given a
% suitable quadrature rule.
qr = QuadratureRule.ofOrder(1, '2D');
val = integrateElement(f, qr);
fprintf('Approximate integral of f: %f\n', sum(val))

% There is also a function to integrate (normal) jumps over edges.
