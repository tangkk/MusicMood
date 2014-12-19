function output = zeroCrossingRate(input)

% ****************************zeroCrossing Rate**********************************%
        songLength = length(input);
        zeroCrossing = 0;
        for i = 1:1:songLength-1
            zeroCrossing = zeroCrossing + abs(sign(input(i+1)) - sign(input(i)));
        end
        output = zeroCrossing / songLength;