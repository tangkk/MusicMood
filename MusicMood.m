% musicmood.m
% This file extracts several features from a piece of music(audio) and plot
% the features as a polygon.
% The features: [zeroCrossingRate, specCentroid, specRolloff, unitPower, lowEnergyRate, tonalType, barTempo];
% Written by Junqi Deng, July 22nd, 2014
% The University of Hong Kong

% ****************************Audio Read**********************************%
function data = musicmood(path)

info = audioinfo(path);
songLength = info.TotalSamples;
sampleRate = info.SampleRate;
numChannels = info.NumChannels;
section = 5;
close all;

% compute music features every 'section' seconds of a song
zeroCrossing = 0;
power = 0;
spectrumWeight = 0;
spectrumAmpSum = 0;
tonalGravity = [];
tonalType = [];
beatTempo = [];
barTempo = [];

sampleStart = 1;
while sampleStart < songLength
    sampleEnd = min(sampleStart + section*sampleRate - 1, songLength);
    [song,fs] = audioread(path, [sampleStart, sampleEnd]);
    sampleStart = sampleEnd + 1;

    % Start calculating other Features
    sizeSong = size(song);
    if sizeSong(2) == 0
        continue;
    end
    if numChannels == 2
        songMono = (song(:,1) + song(:,2)) / 2;
    else
        songMono = song;
    end
    
    sectionLength = length(songMono);
    
    % ****************************unit power**********************************%
    % note that signal power is proportional to the square of signal amplitude
    power = power + sum(songMono.^2);
    
%     % normalization
%     songMonoMax = max(abs(songMono));
%     songMono = songMono / songMonoMax;

    % play music
    % player = audioplayer(songMono, fs);
    % play(player);
    % stop(player);

    % ****************************zeroCrossing Rate**********************************%
    for i = 1:1:sectionLength-1
        zeroCrossing = zeroCrossing + abs(sign(songMono(i+1)) - sign(songMono(i)));
    end

    % ****************************fft transform**********************************%
    NFFT = 2^nextpow2(sectionLength);
    fftSong = fft(songMono,NFFT)/sectionLength;
    f = fs/2*linspace(0,1,NFFT/2+1);
    fftAmpSpec = abs(fftSong(1:NFFT/2+1));

%     plot(f,fftAmpSpec)
%     title('Single-Sided Amplitude Spectrum of song(t)')
%     xlabel('Frequency (Hz)')
%     ylabel('|song(f)|')

    % ****************************spectral centroid**********************************%
    spectrumWeight = spectrumWeight + f*fftAmpSpec;
    spectrumAmpSum = spectrumAmpSum + sum(fftAmpSpec);
    
    % ****************************tonality**********************************%
    tonal = tonality(f, fftAmpSpec);
    tonalType = [tonalType tonal(1)];
    tonalGravity = [tonalGravity tonal(2)];
    
    % ****************************tempo**********************************%
    tmp = tempo(song, fs);
    beatTempo = [beatTempo tmp(1)];
    barTempo = [barTempo tmp(2)];
end
zeroCrossingRate = zeroCrossing / songLength;
unitPower = power / songLength;
specCentroid = spectrumWeight / spectrumAmpSum;
specCentroid = specCentroid / (sampleRate/2); % normalization
tonalT = mode(tonalType);
tonalG = mode(tonalGravity);
beatTempo(beatTempo > 200) = [];
barTempo(barTempo > 50) = [];
bpm = mode(beatTempo);
barpm = mode(barTempo);

%**********************Low Energy****************************%
lowEnergy = 0;
sampleStart = 1;
while sampleStart < songLength
    sampleEnd = min(sampleStart + section*sampleRate - 1, songLength);
    [song,fs] = audioread(path, [sampleStart, sampleEnd]);
    sampleStart = sampleEnd + 1;

    % Start calculating other Features
    sizeSong = size(song);
    if sizeSong(2) == 0
        continue;
    end
    if numChannels == 2
        songMono = (song(:,1) + song(:,2)) / 2;
    else
        songMono = song;
    end
    
    sectionLength = length(songMono);
    % ****************************low energy rate**********************************%
    for i = 1:1:sectionLength
        if songMono(i)^2 < unitPower
            lowEnergy = lowEnergy + 1;
        end
    end
    
end
lowEnergyRate = lowEnergy / songLength;

data = [zeroCrossingRate, specCentroid, unitPower, lowEnergyRate, tonalT, bpm];
display(zeroCrossingRate);
display(specCentroid);
display(unitPower);
display(lowEnergyRate);
display(pitch2name(tonalT));
display(pitch2name(tonalG));
display(bpm);

% ****************************plot**********************************%
% figure;
% bar(chroma);
% set(gca, 'XTickLabel', {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#','B'});
% 
% figure;
% bar(tonalMajorScale);
% set(gca, 'XTickLabel', {'CM', 'C#M', 'DM', 'D#M', 'EM', 'FM', 'F#M', 'GM', 'G#M', 'AM', 'A#M','BM'});
% 
% figure;
% bar(tonalMinorScale);
% set(gca, 'XTickLabel', {'Cm', 'C#m', 'Dm', 'D#m', 'Em', 'Fm', 'F#m', 'Gm', 'G#m', 'Am', 'A#m','Bm'});

% plot bar graph function
% figure;
% bar(data);
% set(gca, 'XTickLabel', {'zeroCrossingRate', 'specCentroid', 'unitPower', 'lowEnergyRate', 'tonal(1)', 'normalized bpm'});

% plot polar graph function
% for i=1:1:length(data)
%     theta(i) = (2*pi/length(data))*i;
% end
% figure;
% polar(theta, data);

% close all;
% fileID = fopen('a.txt','w');
% fprintf(fileID, '%f,',songAutoCorrSec);
