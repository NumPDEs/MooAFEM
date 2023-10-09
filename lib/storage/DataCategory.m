classdef DataCategory
    enumeration
        ERROR(@loglog, 'Convergence history plot', 'error',  'northeast')
        ABSOLUTE(@semilogx, 'Value plot', 'value',  'northeast')
        TIME(@loglog, 'Time plot', 'runtime',  'northwest')
    end

    properties (SetAccess=immutable)
        plotFunction function_handle
        title string
        yLabel string
        legendLocation string
    end

    methods
        function obj = DataCategory(plotFunction, title, yLabel, legendLocation)
            obj.plotFunction = plotFunction;
            obj.title = title;
            obj.yLabel = yLabel;
            obj.legendLocation = legendLocation;
        end
    end
end