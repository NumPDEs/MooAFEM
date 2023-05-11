function retval = isOctave()
%%ISOCTAVE checks whether this function is called from Octave
%
%   Source:
%
%   https://octave.org/doc/v5.1.0/How-to-Distinguish-Between-Octave-and-Matlab.html

    persistent cacheval;  % speeds up repeated calls

    if isempty (cacheval)
        cacheval = (exist ("OCTAVE_VERSION", "builtin") > 0);
    end

    retval = cacheval;
end
