function writevideoMP4(u,name,fps)        

video_3D = VideoWriter(name,'Uncompressed AVI');
video_3D.FrameRate = fps;

open(video_3D);

for tt=1:size(u,4)-1
    %h = figure(1);
    %imshow(u(:,:,:,tt),[0,1],'border','tight');
    %axis off
    %truesize
    writeVideo(video_3D,u(:,:,:,tt)); 
end
close(video_3D);

