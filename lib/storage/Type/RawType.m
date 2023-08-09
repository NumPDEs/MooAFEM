classdef RawType
    enumeration
        INT('d')
        FLOAT('.5e')
        TEXT('s')
    end

    properties (SetAccess=immutable, GetAccess=public)
        formatSpec
        
    end

    methods
        function obj = RawType(formatSpec)
            obj.formatSpec = formatSpec;
        end
    end
end