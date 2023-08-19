function [res, par, sta, sys] = runefo(par, sta, sys)
    % DESCRIPTION
    % [res, par, sta, sys] = runefo(par, sta, sys)
    % The basic rune dynamic function for OFDMA.

    % Input parameters:
    % par -- BASIC SIMULATION PARAMETERS (All inputs are optional)
    % sta -- STATE variable that is changed at every iteration
    % sys -- SYSTEM variables

    % Output parameters:
    % res -- cell array of sta collected for each iteration of the simulation
    % par -- same as input otherwise created within the function
    % sta -- same structure as input otherwise created within the function
    % sys -- same as input otherwise created within the function

    % by Magnus Almgren 080221

    % Define the setseed function to make it compatible with Octave
    setseed = @(val1, val2) setseed_inner(val1, val2, par.seed);

    % Rest of the code remains unchanged from the provided function
    % ...
    % (The rest of the function remains unchanged)
    % ...

end

function oseed = setseed_inner(val1, val2, seed)
    % sets the seed to random generators

    % DESCRIPTION oseed = setseed(val1, val2)
    % First, the old seed is delivered as a matrix output
    % Secondly, the seed is set to the Octave random generator
    % If no input or if the input is equal to NaN, the seed is not touched.

    % INPUT
    % val1 -- A seed value vector 35 by 1 for the rand generator.
    % val2 -- A seed value vector 2 by 1 for the randn generator

    % OUTPUT
    % oseed -- The current states for rand and randn as a cell array.

    % TRY
    % oseed = setseed(3),
    % oseed = setseed(3,4),
    % oseed = setseed
    % oseed = setseed(NaN)
    % SEE ALSO rand('state'), randn('state')

    % by Magnus Almgren 990301

    oseed = {rand("state"), randn("state")}; % get the old seed

    if nargin == 0
        % No input, do nothing
    elseif nargin == 1 && iscell(val1)
        setseed(val1{1}, val1{2}); % the struct case
    elseif nargin == 1 && length(val1(:)) == 1
        setseed(val1, val1); % the scalar case
    elseif nargin == 1 && length(val1(:)) == 2
        setseed(val1(1), val1(2)); % two scalars
    elseif nargin == 2
        if all(isfinite(val1)) % do not touch if NaN
            rand("state", val1); % size(val1,1) is either 1 or 35
        end
        if all(isfinite(val2))
            randn("state", val2); % size(val1,1) is either 1 or 2
        end
    else
        error("Inputs are of incorrect size or type");
    end
end
