function plotTriangulationToFile(obj, jLevel)
%%PLOTTRIANGULATIONTOFILE creates (a series of) plots for triangulations
%and stores it to a file, filename is generated automatically from
%information in LevelData object
%   ax = PLOTTRIANGULATIONTOFILE(obj) plots all levels
%   ax = PLOTTRIANGULATIONTOFILE(obj, ':')
%   ax = PLOTTRIANGULATIONTOFILE(obj, jLevel)
%   ax = PLOTTRIANGULATIONTOFILE(obj, 'end')
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


    % Proceed input
    if nargin < 2
        jLevel = 1:obj.nLevel;
    elseif strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    elseif strcmp(jLevel, 'end')
        jLevel = obj.nLevel;
    end

    % Create problem- and method-specific folder
    ensureFolderExists(obj.foldername);

    % Create new figure
    h = createStandardFigure();

    if length(jLevel) == 1
        % Create single image
        h = createStandardisedFigure();
        obj.plotLevelTriangulation(gca(), jLevel);
        % Export plot
        print(h, '-dpng', '-r600', [obj.foldername, '/', obj.filename, '.png']);
    else
        % Set duration of video in seconds
        DURATION = 20;
        % Create video as series of images
        if isOctave()
            createTriangulationGIF(obj, h, jLevel, DURATION);
        else
            createTriangulationAVI(obj, h, jLevel, DURATION);
        end
    end

    % Close figure
    close(h);
end


function createTriangulationAVI(obj, h, jLevel, duration)
%%CREATETRIANGULATIONAVI creates a video of the mesh refinement in AVI format

    % Set video quality parameter
    FPS = 30;

    % Create video writer
    videowriter = VideoWriter([obj.foldername, '/', obj.filename, '.avi']);
    set(videowriter, 'FrameRate', FPS);
    open(videowriter);

    % Compute time per frame to obtain a fixed overall duration
    framesPerLevel = fix(FPS * duration / length(jLevel));
    ax = gca();
    set(ax, 'NextPlot', 'replaceChildren');

    % Plot every frame
    for k = 1:length(jLevel)
        obj.plotLevelTriangulation(ax, jLevel(k));
        for l = 1:framesPerLevel
            writeVideo(videowriter, getframe(h));
        end
    end

    % Close video writer
    close(videowriter);
end


function createTriangulationGIF(obj, h, jLevel, duration)
%%CREATETRIANGULATIONGIF creates a video of the mesh refinement in GIF format

    error('not working yet')

    ax = gca();

    % Determine time for delay in GIF
    delayTime = duration / length(jLevel);

    % Define path
    filepath = [obj.foldername, '/', obj.filename, '.gif'];

    % Assign plot to a frame
    frame = getframe(h);

    % Convert frame to RGB image (3 dimensional)
    im = frame2im(frame);

    % Transform RGB samples to 1 dimension with a color map "cm".
    [imind, cm] = rgb2ind(im);

    % Compute time per frame to obtain a fixed overall duration
    imwrite(imind, cm, filepath, 'gif', ...
            'DelayTime', delayTime, 'Compression' , 'lzw');

    % Plot every frame
    for k = 1:length(jLevel)
        obj.plotLevelTriangulation(ax, jLevel(k));
            % Assign plot to a frame
    frame = getframe(h);

    % Convert frame to RGB image (3 dimensional)
    im = frame2im(frame);

    % Transform RGB samples to 1 dimension with a color map "cm".
    [imind, cm] = rgb2ind(im);

        % Create GIF file
        % Add each new plot to GIF
        imwrite(imind, cm, filepath, 'gif', ...
                'WriteMode', 'append', 'DelayTime', delayTime, ...
                'Compression', 'lzw');
    end
end
