% MusicMood.m
% This file extracts several features from a piece of music(audio) and plot
% the features as a polygon.
% The features: [zeroCrossingRate, specCentroid, specRolloff, unitPower, lowEnergyRate, tonalType, barTempo];
% Written by Junqi Deng, July 22nd, 2014
% The University of Hong Kong

% ****************************Audio Read**********************************%
 [song,fs] = audioread('haoting.mp3');

songLeft = song(:,1);
songRight = song(:,2);
songMono = (songLeft + songRight) / 2;
songLength = length(songMono);
songSub = songLeft - songRight;

% normalization
songMonoMax = max(abs(songMono));
songMono = songMono / songMonoMax;
songSubMax = max(abs(songSub));
songSub = songSub / songSubMax;

% play music
player = audioplayer(songMono, fs);
playerSub = audioplayer(songSub, fs);
% play(player);
% stop(player);

% ****************************zeroCrossing Rate**********************************%
zeroCrossing = 0;
for i = [1:1:songLength-1]
    zeroCrossing = zeroCrossing + abs(sign(songMono(i+1)) - sign(songMono(i)));
end

zeroCrossingRate = zeroCrossing / songLength;

% ****************************unit power**********************************%
% note that signal power is proportional to the square of signal amplitude
rmsPower = rms(songMono);
unitPower = rmsPower.^2;

% ****************************low energy rate**********************************%
lowEnergy = 0;
for i = [1:1:songLength]
    if songMono(i).^2 < unitPower
        lowEnergy = lowEnergy + 1;
    end
end
lowEnergyRate = lowEnergy ./ songLength;

% ****************************tempo**********************************%
% use auto-correlation method
songDif = songRight - songLeft; %vocal removal
songDifMax = max(songDif);
songDif = songDif / songDifMax;
songDifAbs = abs(songDif);
downSampleRate = 32;
songDifAbs = downsample(songDifAbs, downSampleRate);
songDif = downsample(songDif, downSampleRate);
%songAutoCorr = autocorr(songDif, length(songDif)-1);
songAutoCorr = autocorr(songDifAbs, length(songDifAbs)-1);
songAutoCorr(songAutoCorr < 0) = 0;
beginShift = round(fs/downSampleRate/4);
fourSeconds = round(4*fs/downSampleRate);
songAutoCorrSec = songAutoCorr(beginShift:fourSeconds);
songAutoCorrSec = medfilt1(songAutoCorrSec, 50);

p = polyfit((1:numel(songAutoCorrSec))', songAutoCorrSec, 5);
pval = polyval(p,(1:numel(songAutoCorrSec))');
songAutoCorrSec = songAutoCorrSec - pval;
songAutoCorrSec = songAutoCorrSec / max(songAutoCorrSec);

minPeak = 0.5;
for i=[1:1:5]
    [pks,locs] = findpeaks(songAutoCorrSec,'MINPEAKHEIGHT',minPeak, 'MINPEAKDISTANCE', beginShift,'SORTSTR', 'descend');
    if length(locs) < 3
        minPeak = minPeak - 0.1;
    end
end
locs = locs + beginShift;

figure;
bpm = 60./([1:1:length(songAutoCorr)]./(fs/downSampleRate));
plot(bpm(beginShift:fourSeconds),songAutoCorrSec);
%plot((fs/downSampleRate/3:2*fs),songAutoCorr(fs/downSampleRate/3:2*fs));
title('sectionAutoCorr(t) against BPM');
xlabel('BPM');
ylabel('sectionAutoCorr(t)');

bpmPeaks = bpm(locs);
bpmPeaks = round(bpmPeaks);
bpmPeaks = sort(bpmPeaks, 'descend');
% maybe only suitable for 4/4 time songs
if length(bpmPeaks)>=3
    barTempo = bpmPeaks(length(bpmPeaks));
end
if length(bpmPeaks)>=1
    beatTempo = bpmPeaks(1);
end
% assume the tempo does not change
% songDif = songRight - songLeft; %vocal removal
% songDifMax = max(songDif);
% songDif = songDif / songDifMax;
% songDifAbs = abs(songDif);
% songDifAbs = abs(songDif);
% downSampleRate = 1;
% songDifAbs = downsample(songDifAbs, downSampleRate);
% songDif = downsample(songDif, downSampleRate);
% 
% % use cross-correlation method
% songDifSection = songDif(round(length(songDif)/4) : round(length(songDif)/4) + round(4*fs/downSampleRate));
% songDifRemain = songDif(round(length(songDif)/4) + round(4*fs/downSampleRate) : length(songDif));
% sectionXCorr = xcorr(songDifSection,songDifRemain, round(20*fs/downSampleRate));
% 
% songDifAbsSection = songDifAbs(round(length(songDifAbs)/4) : round(length(songDifAbs)/4) + round(4*fs/downSampleRate));
% songDifAbsRemain = songDifAbs(round(length(songDifAbs)/4) + round(4*fs/downSampleRate) : length(songDifAbs));
% sectionAbsXCorr = xcorr(songDifAbsSection,songDifAbsRemain, round(20*fs/downSampleRate));
% 
% filterCoff = fir1(50,.4);     % 50th-order linear-phase FIR filter.
% hd = dfilt.dffir(filterCoff);    % Direct-form FIR implementation.
% sectionXCorr = filter(hd,sectionXCorr);
% 
% filterCoff = fir1(50,.4);     % 50th-order linear-phase FIR filter.
% hd = dfilt.dffir(filterCoff);    % Direct-form FIR implementation.
% sectionAbsXCorr = filter(hd,sectionAbsXCorr);

% figure;
% NFFT = 2^nextpow2(length(sectionXCorr));
% fftx = fft(sectionXCorr,NFFT)/length(sectionXCorr);
% f = fs/2/downSampleRate*linspace(0,1,NFFT/2+1);
% plot(f,abs(fftx(1:NFFT/2+1)));
% title('Single-Sided Amplitude Spectrum of sectionXCorr(t)');
% xlabel('Frequency (Hz)');
% ylabel('|sectionXCorr(f)|');
% 
% figure;
% NFFT = 2^nextpow2(length(sectionAbsXCorr));
% fftx = fft(sectionAbsXCorr,NFFT)/length(sectionAbsXCorr);
% f = fs/2/downSampleRate*linspace(0,1,NFFT/2+1);
% plot(f,abs(fftx(1:NFFT/2+1)));
% title('Single-Sided Amplitude Spectrum of sectionXCorr(t)');
% xlabel('Frequency (Hz)');
% ylabel('|sectionXCorr(f)|');

% % use multi-fold autocorrelation
% filterCoff = fir1(50,.4);     % 50th-order linear-phase FIR filter.
% hd = dfilt.dffir(filterCoff);    % Direct-form FIR implementation.
% songAutoCorr = filter(hd,songAutoCorr);
% songAutoAutoCorr = autocorr(songAutoCorr, length(songAutoCorr)-1);
% songAutoAutoAutoCorr = autocorr(songAutoAutoCorr, length(songAutoAutoCorr)-1);
% songAutoAutoAutoAutoCorr = autocorr(songAutoAutoAutoCorr, length(songAutoAutoAutoCorr)-1);
% 
% % fft of the above multi-auto-corr
% figure;
% NFFT = 2^nextpow2(length(songAutoAutoAutoAutoCorr));
% fftx = fft(songAutoAutoAutoAutoCorr,NFFT)/length(songAutoAutoAutoAutoCorr);
% f = fs/2/downSampleRate*linspace(0,1,NFFT/2+1);
% plot(f,abs(fftx(1:NFFT/2+1)));
% title('Single-Sided Amplitude Spectrum of songAutoAutoAutoAutoCorr(t)');
% xlabel('Frequency (Hz)');
% ylabel('|songAutoAutoAutoAutoCorr(f)|');
% 
% figure;
% NFFT = 2^nextpow2(length(songAutoCorr));
% fftx = fft(songAutoCorr,NFFT)/length(songAutoCorr);
% f = fs/2/downSampleRate*linspace(0,1,NFFT/2+1);
% plot(f,abs(fftx(1:NFFT/2+1)));
% title('Single-Sided Amplitude Spectrum of songAutoCorr(t)');
% xlabel('Frequency (Hz)');
% ylabel('|songAutoCorr(f)|');


% ****************************fft transform**********************************%
NFFT = 2^nextpow2(songLength);
fftSong = fft(songMono,NFFT)/songLength;
f = fs/2*linspace(0,1,NFFT/2+1);
fftAmpSpec = abs(fftSong(1:NFFT/2+1));

plot(f,fftAmpSpec)
title('Single-Sided Amplitude Spectrum of song(t)')
xlabel('Frequency (Hz)')
ylabel('|song(f)|')

% ****************************spectral centroid**********************************%
spectrumWeight = f*fftAmpSpec;
ampSum = sum(fftAmpSpec);
specCentroid = spectrumWeight ./ ampSum;
specCentroid = specCentroid ./ (fs/2); % normalization

% ****************************spectral roll-off**********************************%
N = length(f);
sumAmp = 0;
for i = [1:1:N]
    sumAmp = sumAmp + fftAmpSpec(i);
    if sumAmp >= 0.85*ampSum
        specRolloff = i*fs/2/N;
        specRolloff = specRolloff ./ (fs/2);
        break;
    end
end

% ****************************tonal gravity and tonal type**********************************%
% the mid feature is chroma of the song
% A = 69, mod (69, 12) = 9, thus
% C 1
% C# 2
% D 3
% D# 4
% E 5
% F 6
% F# 7
% G 8
% G# 9
% A 10
% A# 11
% B 12
for i = [1:1:length(f)]
    if f(i) >= 20
        fm = f(i:length(f));
        fftAmpSpecShort = fftAmpSpec(i:length(fftAmpSpec));
        break;
    end
end
midiNumber = round(12.*log2(fm./440)+69);
midiPitchClass = mod (midiNumber, 12);
midiPitchClass = midiPitchClass + 1;
chroma = [0,0,0,0,0,0,0,0,0,0,0,0];

for i=[1:1:length(midiPitchClass)]
    chroma(midiPitchClass(i)) = chroma(midiPitchClass(i)) + fftAmpSpecShort(i);
end
[chromaMax, tonalMax] = max(chroma);
chroma = chroma/chromaMax;

% C C# D D# E F F# G G# A  A# B
% 1 2  3 4  5 6 7  8 9  10 11 12
% The tonal is calculated according to major and minor scale
tonalMajor = [0,0,0,0,0,0,0,0,0,0,0,0];
tonalMajor(1) = chroma(1)+chroma(5)+chroma(8);
tonalMajor(2) = chroma(2)+chroma(6)+chroma(9);
tonalMajor(3) = chroma(3)+chroma(7)+chroma(10);
tonalMajor(4) = chroma(4)+chroma(8)+chroma(11);
tonalMajor(5) = chroma(5)+chroma(9)+chroma(12);
tonalMajor(6) = chroma(6)+chroma(10)+chroma(1);
tonalMajor(7) = chroma(7)+chroma(11)+chroma(2);
tonalMajor(8) = chroma(8)+chroma(12)+chroma(3);
tonalMajor(9) = chroma(9)+chroma(1)+chroma(4);
tonalMajor(10) = chroma(10)+chroma(2)+chroma(5);
tonalMajor(11) = chroma(11)+chroma(3)+chroma(6);
tonalMajor(12) = chroma(12)+chroma(4)+chroma(7);

tonalMajor7 = [0,0,0,0,0,0,0,0,0,0,0,0];
tonalMajor7(1) = chroma(1)+chroma(5)+chroma(8) + chroma(12);
tonalMajor7(2) = chroma(2)+chroma(6)+chroma(9) + chroma(1);
tonalMajor7(3) = chroma(3)+chroma(7)+chroma(10) + chroma(2);
tonalMajor7(4) = chroma(4)+chroma(8)+chroma(11) + chroma(3);
tonalMajor7(5) = chroma(5)+chroma(9)+chroma(12) + chroma(4);
tonalMajor7(6) = chroma(6)+chroma(10)+chroma(1) + chroma(5);
tonalMajor7(7) = chroma(7)+chroma(11)+chroma(2) + chroma(6);
tonalMajor7(8) = chroma(8)+chroma(12)+chroma(3) + chroma(7);
tonalMajor7(9) = chroma(9)+chroma(1)+chroma(4) + chroma(8);
tonalMajor7(10) = chroma(10)+chroma(2)+chroma(5) + chroma(9);
tonalMajor7(11) = chroma(11)+chroma(3)+chroma(6) + chroma(10);
tonalMajor7(12) = chroma(12)+chroma(4)+chroma(7) + chroma(11);

tonalMinor = [0,0,0,0,0,0,0,0,0,0,0,0];
tonalMinor(1) = chroma(1)+chroma(4)+chroma(8);
tonalMinor(2) = chroma(2)+chroma(5)+chroma(9);
tonalMinor(3) = chroma(3)+chroma(6)+chroma(10);
tonalMinor(4) = chroma(4)+chroma(7)+chroma(11);
tonalMinor(5) = chroma(5)+chroma(8)+chroma(12);
tonalMinor(6) = chroma(6)+chroma(9)+chroma(1);
tonalMinor(7) = chroma(7)+chroma(10)+chroma(2);
tonalMinor(8) = chroma(8)+chroma(11)+chroma(3);
tonalMinor(9) = chroma(9)+chroma(12)+chroma(4);
tonalMinor(10) = chroma(10)+chroma(1)+chroma(5);
tonalMinor(11) = chroma(11)+chroma(2)+chroma(6);
tonalMinor(12) = chroma(12)+chroma(3)+chroma(7);

tonalMinor7 = [0,0,0,0,0,0,0,0,0,0,0,0];
tonalMinor7(1) = chroma(1)+chroma(4)+chroma(8) + chroma(11);
tonalMinor7(2) = chroma(2)+chroma(5)+chroma(9) + chroma(12);
tonalMinor7(3) = chroma(3)+chroma(6)+chroma(10) + chroma(1);
tonalMinor7(4) = chroma(4)+chroma(7)+chroma(11) + chroma(2);
tonalMinor7(5) = chroma(5)+chroma(8)+chroma(12) + chroma(3);
tonalMinor7(6) = chroma(6)+chroma(9)+chroma(1) + chroma(4);
tonalMinor7(7) = chroma(7)+chroma(10)+chroma(2) + chroma(5);
tonalMinor7(8) = chroma(8)+chroma(11)+chroma(3) + chroma(6);
tonalMinor7(9) = chroma(9)+chroma(12)+chroma(4) + chroma(7);
tonalMinor7(10) = chroma(10)+chroma(1)+chroma(5) + chroma(8);
tonalMinor7(11) = chroma(11)+chroma(2)+chroma(6) + chroma(9);
tonalMinor7(12) = chroma(12)+chroma(3)+chroma(7) + chroma(10);

% C C# D D# E F F# G G# A  A# B
% 1 2  3 4  5 6 7  8 9  10 11 12
% The tonal is calculated according to major and minor scale

tonalMajorScale = [0,0,0,0,0,0,0,0,0,0,0,0];
tonalMajorScale(1) = chroma(1)+chroma(3)+chroma(5)+chroma(6)+chroma(8)+chroma(10)+chroma(12);
tonalMajorScale(2) = chroma(2)+chroma(4)+chroma(6)+chroma(7)+chroma(9)+chroma(11)+chroma(1);
tonalMajorScale(3) = chroma(3)+chroma(5)+chroma(7)+chroma(8)+chroma(10)+chroma(12)+chroma(2);
tonalMajorScale(4) = chroma(4)+chroma(6)+chroma(8)+chroma(9)+chroma(11)+chroma(1)+chroma(3);
tonalMajorScale(5) = chroma(5)+chroma(7)+chroma(9)+chroma(10)+chroma(12)+chroma(2)+chroma(4);
tonalMajorScale(6) = chroma(6)+chroma(8)+chroma(10)+chroma(11)+chroma(1)+chroma(3)+chroma(5);
tonalMajorScale(7) = chroma(7)+chroma(9)+chroma(11)+chroma(12)+chroma(2)+chroma(4)+chroma(6);
tonalMajorScale(8) = chroma(8)+chroma(10)+chroma(12)+chroma(1)+chroma(3)+chroma(5)+chroma(7);
tonalMajorScale(9) = chroma(9)+chroma(11)+chroma(1)+chroma(2)+chroma(4)+chroma(6)+chroma(8);
tonalMajorScale(10) = chroma(10)+chroma(12)+chroma(2)+chroma(3)+chroma(5)+chroma(7)+chroma(9);
tonalMajorScale(11) = chroma(11)+chroma(1)+chroma(3)+chroma(4)+chroma(6)+chroma(8)+chroma(10);
tonalMajorScale(12) = chroma(12)+chroma(2)+chroma(4)+chroma(5)+chroma(7)+chroma(9)+chroma(11);

tonalMinorScale = [0,0,0,0,0,0,0,0,0,0,0,0];
tonalMinorScale(1) = chroma(1)+chroma(3)+chroma(4)+chroma(6)+chroma(8)+chroma(9)+chroma(11);
tonalMinorScale(2) = chroma(2)+chroma(4)+chroma(5)+chroma(7)+chroma(9)+chroma(10)+chroma(12);
tonalMinorScale(3) = chroma(3)+chroma(5)+chroma(6)+chroma(8)+chroma(10)+chroma(11)+chroma(1);
tonalMinorScale(4) = chroma(4)+chroma(6)+chroma(7)+chroma(9)+chroma(11)+chroma(12)+chroma(2);
tonalMinorScale(5) = chroma(5)+chroma(7)+chroma(8)+chroma(10)+chroma(12)+chroma(1)+chroma(3);
tonalMinorScale(6) = chroma(6)+chroma(8)+chroma(9)+chroma(11)+chroma(1)+chroma(2)+chroma(4);
tonalMinorScale(7) = chroma(7)+chroma(9)+chroma(10)+chroma(12)+chroma(2)+chroma(3)+chroma(5);
tonalMinorScale(8) = chroma(8)+chroma(10)+chroma(11)+chroma(1)+chroma(3)+chroma(4)+chroma(6);
tonalMinorScale(9) = chroma(9)+chroma(11)+chroma(12)+chroma(2)+chroma(4)+chroma(5)+chroma(7);
tonalMinorScale(10) = chroma(10)+chroma(12)+chroma(1)+chroma(3)+chroma(5)+chroma(6)+chroma(8);
tonalMinorScale(11) = chroma(11)+chroma(1)+chroma(2)+chroma(4)+chroma(6)+chroma(7)+chroma(9);
tonalMinorScale(12) = chroma(12)+chroma(2)+chroma(3)+chroma(5)+chroma(7)+chroma(8)+chroma(10);

% tonalMajorSum = (tonalMajor + tonalMajor7)/2 + tonalMajorScale;
% tonalMinorSum = (tonalMinor + tonalMinor7)/2 + tonalMinorScale;
tonalMajorSum = tonalMajorScale;
tonalMinorSum = tonalMinorScale;
[maxMajor, majorIdx] = max(tonalMajorSum);
[maxMinor, minorIdx] = max(tonalMinorSum);
tonal = [0,0];
% if maxMajor == maxMinor
    if chroma(majorIdx) == chroma(minorIdx)
        tonal = [0,0]; % none key
    else
        if chroma(majorIdx) > chroma(minorIdx)
            tonal = [1, majorIdx]; % major key
        else
            tonal = [-1, minorIdx]; % minor key
        end
    end
% else
%     if maxMajor > maxMinor
%         tonal = [1, majorIdx]; % major key
%     else
%         tonal = [-1, minorIdx]; % minor key
%     end
% end

% ****************************plot**********************************%
figure;
bar(chroma);
set(gca, 'XTickLabel', {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#','B'});

figure;
bar(tonalMajorScale);
set(gca, 'XTickLabel', {'CM', 'C#M', 'DM', 'D#M', 'EM', 'FM', 'F#M', 'GM', 'G#M', 'AM', 'A#M','BM'});

figure;
bar(tonalMinorScale);
set(gca, 'XTickLabel', {'Cm', 'C#m', 'Dm', 'D#m', 'Em', 'Fm', 'F#m', 'Gm', 'G#m', 'Am', 'A#m','Bm'});

% plot bar graph function
data = [zeroCrossingRate, specCentroid, unitPower, lowEnergyRate, tonal(1), barTempo/200];
figure;
bar(data);
set(gca, 'XTickLabel', {'zeroCrossingRate', 'specCentroid', 'unitPower', 'lowEnergyRate', 'tonal(1)', 'normalized bpm'});

% plot polar graph function
for i=[1:1:length(data)]
    theta(i) = (2*pi/length(data))*i;
end
figure;
polar(theta, data);

% close all;
% fileID = fopen('a.txt','w');
% fprintf(fileID, '%f,',songAutoCorrSec);
