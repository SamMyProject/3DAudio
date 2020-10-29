clear;
close all;
format long
addpath(genpath('/home/mvnl/Downloads/API_MO-master'));
%%
SOFAstart;

% Subject index of the file to convert
if ~exist('subject','var'), subject=80; end;
% File name of the RIEC file
RIECfn='/home/mvnl/Downloads/RIEC_hrtf_all/RIEC_hrir_subject_';
RIECfn='/home/mvnl/Downloads/qu_kemar_anechoic_3m.sofa';
% SourcePosition
azim=20; %azimuth
elev=0;  %elevation

%% Path definitions
SOFAfile=[RIECfn sprintf('%03d',subject) '.sofa'];
SOFAfile=RIECfn;
%% Loading the full object
disp(['Loading full object: ' SOFAfile]);
Obj=SOFAload(SOFAfile);

%% Get index of measurements with the same directions
% idx=find(Obj.SourcePosition(:,1)==azim & Obj.SourcePosition(:,2)==elev);
% 
% idx2=find(Obj.SourcePosition(:,1)==270 & Obj.SourcePosition(:,2)==elev);
% 
% %% Extract and plot the fully loaded data
% RIEC_hrirL=squeeze(Obj.Data.IR(idx,1,:)); %left ear
% RIEC_hrirR=squeeze(Obj.Data.IR(idx,2,:)); %right ear

%%%
[y1,Fs1] = audioread('music/There For You.mp3');
audiowrite('music/result.wav',y1,44100,'BitsPerSample',32);

[y,Fs] = audioread('music/result.wav');
% temp=y;
% y = lowpass(y,1600,Fs);
% temp=temp-y;
% temp=temp';

% y = resample(y,48000,44100);
% Fs=48000;
% audiowrite('left.wav',y(:,1),Fs);
% audiowrite('right.wav',y(:,2),Fs);

lc=[];
rc=[];
d1=1;
d2=1;
n=1;
c=1;
re=511;
azim=0;
for i=3:length(y(:,1)')/(Fs)
        idx=find(Obj.SourcePosition(:,1)==mod(azim,360) & Obj.SourcePosition(:,2)==elev);
        idx2=find(Obj.SourcePosition(:,1)==mod(azim,360));
        idx2=mod(mod(azim,360)+181,360);
        l=squeeze(Obj.Data.IR(idx2,1,:)); %left ear
        r=squeeze(Obj.Data.IR(idx2,2,:)); %right ear
        if (d1==1) 
            azim=azim+15;
        elseif(d1==0)
            azim=azim-15;
        end
        
        if (d2==1) 
            c=c+1;
        elseif(d2==0)
            c=c-1;
        end
        e=i*(Fs);
        if e>length(y(:,1)') 
            break;
        end
        st=e-(Fs*3-1);
%         L=conv(y(st:e,1)',squeeze(l(n,33,:))');
%         R=conv(y(st:e,2)',squeeze(r(n,33,:))');
%         L=conv(y(st:e,1)',l');
%         R=conv(y(st:e,2)',r');
%         L=ramp(L(re:Fs+length(l')-re),length(L(re:Fs+length(l')-re)));
%         R=ramp(R(re:Fs+length(r')-re),length(R(re:Fs+length(r')-re)));
        L=ifft(fft(y(st:e,1)').*fft(l',length(y(st:e,1)')));
        R=ifft(fft(y(st:e,2)').*fft(r',length(y(st:e,1)')));
        L=L(44101:88200);
        R=R(44101:88200);

        L=ramp(L,length(L));
        R=ramp(R,length(R));
%         b=find(L==0);
%         L(b)=temp(b)/5;
%         b=find(R==0);
%         R(b)=temp(b)/5;
%         L=ramp(L,length(L))+temp(1,st+44100:e-44100)/10;
%         R=ramp(R,length(R))+temp(2,st+44100:e-44100)/10;
          
          
        
        if(azim==360)
            d1=0;
        elseif(azim==0)
            d1=1;
        end
        if(c==25)
            d2=0;
        elseif(c==1)
            d2=1;
        end

        
        lc=[lc L];
        rc=[rc R];
 
end
lc=lc';
rc=rc';
lc = resample(lc,192000,44100);
rc = resample(rc,192000,44100);
lc=lc/max(abs(lc));
rc=rc/max(abs(rc));
a=[lc,rc];

audiowrite('music/result1.wav',a,192000,'BitsPerSample',32);


function sigfix=pha(sig,hrtf)
    h=conv(sig,hrtf);
    org=fft(sig,96000)';
    hfft=fft(h,96000)';
    phase_shift=unwrap(angle(hfft./org));
    sigfix=real(ifft((fft(sig,96000)'.*exp(-1i*phase_shift))'));
end

function m=ramp(sig,len)
    for i=1:100
        sig(i)=sig(i)*(i-1)/100;
        %disp(sig(i));
    end
    
    for i= (len-99):len
        sig(i)=sig(i)*(len-i)/100;
        %disp(sig(i));
    end
    
    m=sig;
    
end
