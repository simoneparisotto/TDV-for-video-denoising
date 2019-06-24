function [uorig,fps,M,N,C,T,I_MAX] = create_video(object)

%% CREATE/LOAD ORIGINAL VIDEO
% uorig  = real data
% u      = noisy data
% params = parameters for TV3D
videomat = ['./dataset/videomat/',object,'.mat'];
if ~exist(videomat,'file')
    [uorig, fps] = create_data(object);
    save(videomat,'uorig','fps');
else
    load(videomat,'uorig','fps');
end
[M,N,C,T] = size(uorig);
I_MAX = 1;

return