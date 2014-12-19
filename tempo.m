function tmp = tempo(song, fs)
 
% ****************************tempo**********************************%
warning('off','all');
% use auto-correlation method
songDif = song(:,1) - song(:,2); %vocal removal
songDifMax = max(songDif);
songDif = songDif / songDifMax;
songDif = abs(songDif);
downSampleRate = 32;
songDif = downsample(songDif, downSampleRate);
songAutoCorr = myAutocorr(songDif, length(songDif)-1);
songAutoCorr(songAutoCorr < 0) = 0;
beginShift = round(fs/downSampleRate/4);
fourSeconds = round(4*fs/downSampleRate);
if length(songAutoCorr) < fourSeconds
    tmp = [0, 0];
    return;
end
songAutoCorrSec = songAutoCorr(beginShift:fourSeconds);
songAutoCorrSec = medfilt1(songAutoCorrSec, 50);
p = polyfit(1:numel(songAutoCorrSec), songAutoCorrSec, 5);
pval = polyval(p,1:numel(songAutoCorrSec));
songAutoCorrSec = songAutoCorrSec - pval;
songAutoCorrSec = songAutoCorrSec / max(songAutoCorrSec);

minPeak = 0.7;
for i=1:1:5
    [pks,locs] = findpeaks(songAutoCorrSec,'MINPEAKHEIGHT',minPeak, 'MINPEAKDISTANCE', beginShift,'SORTSTR', 'descend');
    if length(locs) < 3
        minPeak = minPeak - 0.1;
    end
end
locs = locs + beginShift;

bpm = 60./(linspace(1,length(songAutoCorr),length(songAutoCorr))./(fs/downSampleRate));
% figure;
% songAutoCorrSec(songAutoCorrSec < 0) = 0;
% plot(bpm(beginShift:fourSeconds),songAutoCorrSec);
% title('sectionAutoCorr(t) against BPM');
% xlabel('BPM');
% ylabel('sectionAutoCorr(t)');

bpmPeaks = bpm(locs);
bpmPeaks = round(bpmPeaks);
bpmPeaks = sort(bpmPeaks, 'descend');
% maybe only suitable for 4/4 time songs
if length(bpmPeaks)>=3
    barTempo = bpmPeaks(length(bpmPeaks));
else
    barTempo = 0;
end

if length(bpmPeaks)>=1
    beatTempo = bpmPeaks(1);
else
    beatTempo = 0;
end

tmp = [beatTempo, barTempo];

