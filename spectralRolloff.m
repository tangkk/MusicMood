function specRolloff = spectralRolloff(f, fftAmpSpec)

    % ****************************spectral roll-off**********************************%
    N = length(f);
    sumAmp = 0;
    for i = 1:1:N
        sumAmp = sumAmp + fftAmpSpec(i);
        if sumAmp >= 0.85*ampSum
            specRolloff = i*fs/2/N;
            specRolloff = specRolloff ./ (fs/2);
            break;
        end
    end