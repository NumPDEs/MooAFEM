function propgrp = getPropertyGroups(obj)
    if ~isscalar(obj)
        propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
    else
        % some rudimentary statistics about variable counts
        statList = ["nLevel", "nVariable", "nScalarVariable", "nAbsoluteVariable", "nTimeVariable"];
        statTitle = "Variable Counts";
        statGrp = matlab.mixin.util.PropertyGroup(statList, statTitle);

        % metadata as a struct
        metaDataList = struct();
        for k = obj.metaData.keys'
            metaDataList.(k{1}) = obj.metaData(k{1});
        end
        metaDataTitle = "Metadata";
        metaDataGrp = matlab.mixin.util.PropertyGroup(metaDataList, metaDataTitle);

        % storage related properties
        storeList = ["root", "foldername", "filename"];
        storeTitle = "Storage details";
        storeGrp = matlab.mixin.util.PropertyGroup(storeList, storeTitle);
        propgrp = [statGrp, metaDataGrp, storeGrp];
    end
end