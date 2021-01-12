% addPitch2RhythmStruct

function [RhythmStruct] = addPitch2RhythmStruct(RhythmStruct,data_dir)


rdir=pwd;

for i = 1:numel(RhythmStruct)
    
    t_SAP = (0:length([RhythmStruct(i).rawamplitude])-1)/1000;
    filename = [RhythmStruct(i).wavname];

    cd(data_dir)

    [x,fs] = audioread(filename);

    [f0mtlb,locs] = pitch(x,fs); % TURNS OUT THIS LOOKS PRETTY GOOD. COMPARABLE TO YIN_ESTIMATOR, AND MEANINGFUL. except for length of outputs which is weird. 100Hz only?!.
    %t = (0:length(x)-1)/fs;
    t0mtlb = (locs-1)/fs;

    % make this ms wise like the other vectors. just resample, stupid workaround but whatever.
    P=size(t_SAP,2);
    Q=size(t0mtlb,1);

    pitchvect = resample(f0mtlb,P,Q);
    
    RhythmStruct(i).rawpitch = pitchvect;
    
end

cd(rdir)