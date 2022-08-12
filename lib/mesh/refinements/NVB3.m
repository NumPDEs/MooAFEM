% NVB3 (subclass of NVB) Class realizing the NVB3 refinement rule (basically
%   just an alias for NVB).

classdef NVB3 < NVB   
    %% public methods
    methods (Access=public)
        function obj = NVB3(mesh)
            obj = obj@NVB(mesh);
        end
    end
end
