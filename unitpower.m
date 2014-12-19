function output = unitpower(input)

        % ****************************unit power**********************************%
        % note that signal power is proportional to the square of signal amplitude
        rmsPower = rms(input);
        output = rmsPower.^2;