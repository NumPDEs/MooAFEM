function plotTriangulationToFile(obj, jLevel, opt)
%%PLOTTRIANGULATIONTOFILE creates (a series of) plots for triangulations
%and stores it to a file, filename is generated automatically from
%information in LevelData object
%   ax = PLOTTRIANGULATIONTOFILE(obj) plots all levels
%   ax = PLOTTRIANGULATIONTOFILE(obj, jLevel) plots levels specified in jLevel
%   ax = PLOTTRIANGULATIONTOFILE(obj, jLevel, opt) additionally specifies options:
%      - timePerLevel (default 1s)
%
%   See also LevelData/plotTriangulation

% Copyright 2023 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

    arguments
        obj
        jLevel (1,:) {mustBeIndexVector} = ':'
        opt.timePerLevel {mustBeNumeric} = 1
    end

    if strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    end

    ensureFolderExists(obj.foldername);
    h = createStandardFigure();

    if length(jLevel) == 1
        % Create single image
        coordinates = obj.level(jLevel).coordinates;
        elements = obj.level(jLevel).elements;
        mesh = Mesh(coordinates, elements, {});
        plot(mesh);
        print(h, '-dpng', '-r600', obj.foldername + '/' + obj.filename + '.png');
    else
        % Create video as series of images
        createTriangulationAVI(obj, h, jLevel, opt.timePerLevel);
    end

    close(h);
end


function createTriangulationAVI(obj, h, jLevel, secondsPerFrame)
    % video quality parameter
    FPS = 1 / secondsPerFrame;

    % Create video writer
    videowriter = VideoWriter(obj.foldername + '/' + obj.filename + '.avi');
    set(videowriter, 'FrameRate', FPS);
    open(videowriter);

    ax = gca();
    set(ax, 'NextPlot', 'replaceChildren');

    % Plot every frame
    for k = jLevel
        coordinates = obj.level(k).coordinates;
        elements = obj.level(k).elements;
        mesh = Mesh(coordinates, elements, {});
        plot(mesh, 'FaceAlpha', 0.0);
        writeVideo(videowriter, getframe(h));
    end

    close(videowriter);
end
