function [uorig, fps] = create_data(object)

%% LOAD VIDEO
switch object
    
    case {'function','function_colour'}
        M = 120;
        N = 120;
        if strcmp(object,'function')
            C = 1;
        else
            C = 3;
        end
        T = 120;
        x = (0:1:(N-1))/N;
        y = ((M-1):-1:0)/M;
        [X,Y] = meshgrid(x,y);
        uorig =  zeros(M,N,C,T);
        temp  =  zeros(M,N,1,T);
        
        c = parula(120);
        
        for t=1:T
            temp(:,:,1,t) = franke(X+0.1*cos(2*pi*t/T),Y+0.3*sin(2*pi*t/T));
        end
        
        temp = (temp-min(temp(:)))./(max(temp(:))-min(temp(:))); % rescale in [0,1]
        
        
        if strcmp(object,'function')
           uorig = temp;
        else
            for t=1:T
                uorig(:,:,1,t) = c(t,1)*temp(:,:,1,t);
                uorig(:,:,2,t) = c(t,2)*temp(:,:,1,t);
                uorig(:,:,3,t) = c(t,3)*temp(:,:,1,t);
            end
        end
        
        [M,N,C,T] = size(uorig);
         
           
    case {'xylophone','xylophone_colour'}
        v = VideoReader('xylophone.mp4');
        M = v.Height;
        N = v.Width;
        C = 3;
        T = uint64(v.FrameRate*v.Duration);
        rescaling = 1;
        uorig = imresize(zeros(M,N,C,T,'uint8'),rescaling);
        
        t = 1;
        while hasFrame(v)
            temp = imresize(readFrame(v),rescaling);
            uorig(:,:,:,t) = temp;
            t=t+1;
        end
        %T = 120;
        T = size(uorig,4);
        uorig = double(uorig(:,:,:,1:T))/255; % select only first T frames of the video
        if strcmp(object,'xylophone')
            uorig = 0.2989*uorig(:,:,1,:) + 0.5870*uorig(:,:,2,:) + 0.1140*uorig(:,:,3,:);
        end
        [M,N,C,T] = size(uorig);
        
        
    case {'drop','drop_colour'} 
        v = VideoReader('./dataset/drop.mov');
        M = v.Height;
        N = v.Width;
        C = 3;
        T = uint64(v.FrameRate*v.Duration);
        rescaling = 0.6;
        uorig = imresize(zeros(M,N,C,T,'uint8'),rescaling);
        
        t = 1;
        while hasFrame(v)
            temp = imresize(readFrame(v),rescaling);
            uorig(:,:,:,t) = temp;
            t=t+1;
        end
        T = 120;
        %T = 40;
        uorig = double(uorig(100:279,255:574,:,40:40+T-1))/255; % select only first T frames of the video
        if strcmp(object,'drop')
            uorig = 0.2989*uorig(:,:,1,:) + 0.5870*uorig(:,:,2,:) + 0.1140*uorig(:,:,3,:);
        end
        [M,N,C,T] = size(uorig);
        
        
    case {'flower','bus','salesman','coastguard','missa','bicycle','foreman','tennis'}
        v = VideoReader(['g',object,'.avi']);
        M = v.Height;
        N = v.Width;
        C = 1;
        T = uint64(v.FrameRate*v.Duration);
        
        rescaling = 1;
        uorig = imresize(zeros(M,N,C,T),rescaling);
        
        t=1;
        while hasFrame(v)
            temp = imresize(readFrame(v),rescaling);
            uorig(:,:,:,t) = temp;
            t=t+1;
        end
        %T = 120;
        T = size(uorig,4);
        uorig = double(uorig(:,:,:,1:T))/255; % select only first T frames of the video
        [M,N,C,T] = size(uorig);
        
       
end

if exist('v','var')
    fps   = v.FrameRate;
else
    fps   = 5;
end