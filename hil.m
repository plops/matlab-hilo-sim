% prepare the test object

%%
NA=1.4;
n=128;
nmperpixel=100;
sz=2;
znmperpixel=sz*nmperpixel;
%% vector psf for excitation illumination
asf=squeeze(kSimPSF({'lambdaEx',473;'na',NA;'ri',1.518;...
    'sX',n;'sY',n;'sZ',sz*n;...
    'scaleX',nmperpixel;'scaleY',nmperpixel;'scaleZ',znmperpixel;...
    'computeASF',1}));


%% grating of period P
% one pixel is 100nm so the period is P*100nm
P=6;
grat2d=(mod(xx(n,n),P)>(floor(P/2)-1))*1000;
%% fill a 3d grating
grat=newim(n,n,sz*n);
grat(:,:,(sz*n)/2)=grat2d(:,:);
kgrat=ft(grat);

%% coherent imaging
imgratx=ift(kgrat.*ft(squeeze(asf(:,:,:,0))));
imgraty=ift(kgrat.*ft(squeeze(asf(:,:,:,1))));
imgratz=ift(kgrat.*ft(squeeze(asf(:,:,:,2))));
imgrat=abs(imgratx).^2+abs(imgraty).^2+abs(imgratz).^2

%% write a x-z section into grating_xz.eps
imgrat_scaled=affine_trans(squeeze(imgrat(:,64,:)),[1,.5],[0,0],0);
figure1 = figure('Visible','off');
axes1 = axes('Parent',figure1,...
    'YTickLabel',{'-6.4','-3.2','0','3.2','6.4'},...
    'YTick',[0.5 32.5 64.5 96.5 128.5],...
    'YDir','reverse',...
    'XTickLabel',{'2','4','6','8','10','12'},...
    'TickDir','out',...
    'Layer','top',...
    'DataAspectRatio',[1 1 1],...
    'CLim',[0 1]);
xlim(axes1,[0.5 128.5]);
ylim(axes1,[0.5 128.5]);
xlabel(axes1,'x/\mum');
ylabel(axes1,'z/\mum');
box(axes1,'on');
hold(axes1,'all');
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 7, 7]);

image(double(imgrat_scaled(:,64:191)),'Parent',axes1,'CDataMapping','scaled');
colormap('copper');

%print(figure1,'-depsc2','/home/martin/0807/grating_xz.eps');
%% for imaging the fluorescence light
psf=kSimPSF({'lambdaEx',473;'lambdaEm',520;'na',1.4;'ri',1.518;...
        'sX',n;'sY',n;'sZ',sz*n;...
        'scaleX',nmperpixel;'scaleY',nmperpixel;'scaleZ',znmperpixel});
%% ft thereof
kpsf=ft(psf);
%% object plane
%obj=newim(imgrat);
%obj(13:114,13:93,(sz*n)/2)=4;
%% as an object I want a hollow sphere
% I define it in k-space
kobj=sinc(rr(kpsf)./2).*sinc(rr(kpsf)./0.7);

% this is the sphere
obj=rr(psf).*abs(ift(kobj));
maximum=max(obj);
% in-focus rectangle in right top
obj(83:114,23:43,(sz*n)/2)=4*maximum;
% in-focus line on the left
obj(21:21,40:90,(sz*n)/2)=12*maximum;

%% shift the object a little bit in z
kobj=ft(obj);
dz=0; % shift in pixels -> 1 equals 100nm
kobj=kobj.*exp(-i*2*pi*zz(kobj,'freq')*dz);
obj=ift(kobj);

%% excited fluorophores
fluo=obj.*imgrat;

%% fluorescence image with structured illumination
strucflimg=ift(ft(fluo).*kpsf);

%% widefield fluorescence image
flimg=ift(ft(obj).*kpsf);

%% extract focal planes
iu=real(squeeze(flimg(:,:,(sz*n)/2)));
in=real(squeeze(strucflimg(:,:,(sz*n)/2)));
%% save defocused structured illumination image
%insb=in*255/max(in);
%insb(sbx:sbx+sbw,sby:sby+sbh)=255;
%writeim(insb,'/home/martin/0807/in100.eps','EPS',0);
%% add scalebar
sbx=100;
sby=110;
sbh=3;
sbw=20; % 2mu
objsb=squeeze(obj(:,:,128));
objsb=objsb*255/max(objsb);
objsb(sbx:sbx+sbw,sby:sby+sbh)=255;
%writeim(objsb,'/home/martin/0807/obj.eps','EPS',0);
insb=in*255/max(in);
insb(sbx:sbx+sbw,sby:sby+sbh)=255;
%writeim(insb,'/home/martin/0807/in.eps','EPS',0);
iusb=iu*255/max(iu);
iusb(sbx:sbx+sbw,sby:sby+sbh)=255;
%writeim(iusb,'/home/martin/0807/iu.eps','EPS',0);
%% project otf along z
skpsf=squeeze(sum(kpsf,[],3));
corr=gaussf((rr(skpsf,'freq')<.42),3)./skpsf; % use this to correct for otf

