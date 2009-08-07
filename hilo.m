%%
n=128;
%% point spread function of our 63x NA1.4 objective blue light in, green light out
nmperpixel=100; % 80
NA=1.4;
psf=kSimPSF({'lambdaEx',473;'lambdaEm',510;'na',NA;'ri',1.518;'sX',n;'sY',n;'sZ',n;'scaleX',nmperpixel;'scaleY',nmperpixel;'scaleZ',2*nmperpixel});
%% ft thereof
kpsf=ft(psf);
%% as an object I want a hollow sphere
% I define it in k-space
kobj=sinc(rr(n,n,n)./2).*sinc(rr(n,n,n)./0.7);

%% this is the sphere
obj=rr(n,n,n).*abs(ift(kobj));
maximum=max(obj);
% in-focus rectangle in right top
obj(83:114,23:43,n/2)=4*maximum;
% in-focus line on the left
obj(21:21,40:90,n/2)=12*maximum;

%% grating of period P
% one pixel is 100nm so the period is 400nm
P=4;
grat2d=(mod(xx(n,n),P)>(floor(P/2)-1))*1000;

%% fill a 3d grating
grat=newim(n,n,n);
grat(:,:,n/2)=grat2d(:,:);
kgrat=ft(grat);

%% coherent illumination
% makes an image of the grating
imgrat=abs(ift(kgrat.*kpsf))^2;

%% shift the object a little bit in z
kobj=ft(obj);
dz=1.6;
kobj=kobj.*exp(-i*2*pi*zz(kobj,'freq')*dz);
obj=ift(kobj);

%% excited fluorophores
fluo=obj.*imgrat;

%% fluorescence image with structured illumination
strucflimg=ift(ft(fluo).*kpsf);

%% widefield fluorescence image
flimg=ift(ft(obj).*kpsf);

%% cut out the plane of the 3d image which is detected by the sensor
s=abs(strucflimg(:,:,n/2));
f=abs(flimg(:,:,n/2));
o=obj(:,:,n/2);
% an in-focus fluorescent plane would give the following image
implane=abs(ift(ft(abs(imgrat)).*kpsf));
il=squeeze(implane(:,:,n/2));
g=squeeze(grat2d(:,:));
%%
asf=squeeze(kSimPSF({'lambdaEx',473;'na',NA;'ri',1.518;'sX',n;'sY',n;'sZ',n;...
    'scaleX',nmperpixel;'scaleY',nmperpixel;'scaleZ',nmperpixel;...
    'computeASF',1}));
%%
psf=sum(abs(asf)^2,[],4);
fasf=dip_fouriertransform(asf,'forward',[1 1 1 0]);
ex=squeeze(s(:,:,:,0));
ey=squeeze(s(:,:,:,1));
ez=squeeze(s(:,:,:,2));
intens=sqrt(ex .* conj(ex) + ey .* conj(ey) + ez .* conj(ez));
%% DampEdge 0.08 2 introduces a column and a row of zeros - not good
iu=real(squeeze(flimg(:,:,n/2)));
in=real(squeeze(strucflimg(:,:,n/2)));