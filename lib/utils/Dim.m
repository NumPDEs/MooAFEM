% Dim (enumeration class) Maps array dimensions to the FEM quantities they
%   represent.

classdef Dim < double
    enumeration
        Vector (1)
        Elements (2)
        QuadratureNodes (3)
    end
end