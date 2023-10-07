% QuadratureRule Contains (barycentric) nodes and weights for numerical
%   quadrature.
%
%   qr = QuadratureRule(bary, weights) constructs QuadratureRule with given
%       nodes (Barycentric) and equally many weights.
%
%   qr = QuadratureRule.ofOrder(p) constructs QuadratureRule with symmetric
%       coordinates up to order 5. For higher orders, non-symmetrical quadrature
%       rules are constructed by the Duffy transform with 1D Gauss-Legendre
%       quadrature rules.
%
%   qr = QuadratureRule.ofOrder(p, '1D') constructs 1D QuadratureRule from
%       Gauss-Legendre quadrature of arbitrary order.

classdef QuadratureRule
    %% properties
    properties (GetAccess=public, SetAccess=protected)
        bary
        weights
        nNodes
        dim
    end
    
    %% methods
    methods (Access=public)
        function obj = QuadratureRule(bary, weights)            
            arguments
                bary (1,1) Barycentric
                weights (1,:) double
            end
            
            obj.bary = bary;
            obj.weights = weights;
            obj.nNodes = bary.nNodes;
            obj.dim = bary.dim;
            
            assertCorrectNumber(bary, weights)
        end
    end
    
    methods(Static)
        % Return quadrature rule of specified order and dimension
        function obj = ofOrder(order, dim)
            arguments
                order {mustBePositive(order)}
                dim {mustBeTextScalar, mustBeMember(dim, {'1D', '2D'})} = '2D'
            end
            
            switch dim
                case '1D'
                    [bary, weights] = quad1d(order);
                    bary = Barycentric1D(bary');
                case '2D'
                    [bary, weights] = quad2d(order);
                    bary = Barycentric2D(bary');
            end
            
            obj = QuadratureRule(bary, weights);
        end

        % Return given quadrature rule or create new one if given rule is
        % an UnspecifiedQR
        function qr = ifEmptyOfOrder(qrOrUnspecified, order, dim)
            arguments
                qrOrUnspecified (1,1) QuadratureRule
                order (1,1) double
                dim {mustBeTextScalar, mustBeMember(dim, {'1D', '2D'})} = '2D'
            end

            if isa(qrOrUnspecified, 'UnspecifiedQR')
                order = max(order, 1);
                qr = QuadratureRule.ofOrder(order, dim);
            else
                qr = qrOrUnspecified;
            end
        end
    end
end

%% local function for correctness assertions
function assertCorrectNumber(bary, weights)
    if bary.nNodes ~= size(weights, Dim.Elements)
        eid = 'QuadratureRule:numberMismatch';
        msg = 'Same number of nodes and weights is required.';
        throwAsCaller(MException(eid,msg))
    end
end

%% local functions generating 1D and 2D quadrature rules
function [bary, weights] = quad1d(order)
    % table of quadrature formulae up to order 9, else compute them on the fly
    npoints = ceil((order+1)/2);
    switch npoints
        case 1
            bary = 0.5;
            weights = 1;
        case 2
            bary =  2.1132486540518712e-1;
            weights = 0.5;
        case 3
            bary = [1.1270166537925831e-1;0.5];
            weights = [5/18;4/9];
        case 4
            bary = [6.9431844202973712e-2; 3.3000947820757187e-1];
            weights = [1.7392742256872693e-1; 3.2607257743127307e-1];
        case 5
            bary = [4.6910077030668004e-2; 2.3076534494715845e-1; 0.5];
            weights = [1.1846344252809454e-1; 2.3931433524968323e-1; 128/450];
        otherwise
            [bary, weights] = lgwt(npoints, 0, 1);
            bary = [flipud(bary),bary];
            return
    end
    
    % generate quadrature points and weights from data above
    weights(npoints:-1:npoints-size(bary,1)+1,1) = weights;
    bary(npoints:-1:npoints-size(bary,1)+1,1) = 1-bary;
    bary = [flipud(bary),bary];
end

function [bary, weights] = quad2d(order)
    % efficient symmetrical quadrature rules for the triangle
    % nodes and weights are originlly taken from
    %
    %   D.A. Dunavant, High degree efficient symmetrical gaussian quadrature
    %   rules for the triangle,  Int J Numer Meth Eng (21) 1985
    %
    %   L. Zhang, T. Cui, and H. Liu, A set of symmetric quadrature rules 
    %   on triangles and tetrahedra, J Comp Math (27) 2009
    %
    % Remark: only quadrature rules with positive weights are listed
    %         the exact formulae of order 5 are also given in
    %         Abramowitz/Stegun: Pocketbook of Mathematical Functions (25.4.63)
    %
    %   This implementation is taken from [Funken, Praetorius, Wissgott; 2011]
    
    if order> 6
        % fallback: Gauss-Legendre + Duffy transform
        [b,w] = quad1d(order+1);
        b = b(:,2);
        weights = 2*asVector(w * (w.*(1-b))');
        bary = [repelem(b, length(b), 1), asVector(b * (1-b)')];
        bary = [1-sum(bary,2),bary];
        return
    end
    
    % table of quadrature formulae on triangle of different order
    wdata = cell(3,1);
    if order <= 1 
        wdata{1} = 1;
    elseif  order <= 2            
        xdata{2} = 2/3;
        wdata{2} = 1/3;
    elseif order <= 4  
        xdata{2}(1:2,1) = 1/9*(1+sqrt(10)-sqrt(38-44*sqrt(2/5))*[1;-1]);
        wdata{2} = (620+sqrt(213125-53320*sqrt(10))*[1;-1])/3720;
    elseif order <= 5 
        wdata{1} = 9/40;
        xdata{2} = (9-2*sqrt(15)*[1;-1])/21;
        wdata{2} = (155+sqrt(15)*[1;-1])/1200;
    elseif order <= 6 
        xdata{2} = [5.0142650965817916e-1; 8.7382197101699554e-1];
        wdata{2} = [1.1678627572637937e-1; 5.0844906370206817e-2];
        xdata{3} = [5.3145049844816947e-2, 3.1035245103378441e-1];
        wdata{3} = [8.2851075618373575e-2] ;
    end
    
    % generate quadrature points and weights from data above
    bary = zeros(0,3); weights=[];
    if ~isempty(wdata{1})
        bary = [1,1,1]/3; 
        weights = wdata{1}; 
    end
    isucc = [2,3,1];
    iprae = [3,1,2];
    if ~isempty(wdata{2})
        points = [xdata{2},(1-xdata{2})/2*[1,1]];
        for i=1:3
            m = [i,isucc(i),iprae(i)];
            bary = [bary;points(:,m)]; 
            weights = [weights;wdata{2}]; 
        end
    end
    if ~isempty(wdata{3})
        points = [xdata{3},1-sum(xdata{3},2)];
        for i=1:3
            m = [i,isucc(i),iprae(i)];
            bary = [bary;points(:,m);points(:,m([1,3,2]))]; 
            weights = [weights;wdata{3};wdata{3}]; 
        end
    end
end
