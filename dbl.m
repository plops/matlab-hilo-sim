run './dbl_g.m';

%% space for the 3d grating
grat=newim(s1,s2,s3);
kgrat=ft(grat);

%% now simulate the coherent imaging into the sample
asf=squeeze(kSimPSF({'lambdaEx',473;'na',1.4;'ri',1.518;...
    'sX',s1;'sY',s2;'sZ',s3;...
    'scaleX',res1;'scaleY',res2;'scaleZ',res3;...
    'computeASF',1}));


%% coherent imaging
kasf0=ft(squeeze(asf(:,:,:,0)));
kasf1=ft(squeeze(asf(:,:,:,1)));
kasf2=ft(squeeze(asf(:,:,:,2)));
imgratx=ift(kgrat.*kasf0);
imgraty=ift(kgrat.*kasf1);
imgratz=ift(kgrat.*kasf2);
imgrat=abs(imgratx).^2+abs(imgraty).^2+abs(imgratz).^2;

%%
% for imaging the fluorescence light
psf=kSimPSF({'lambdaEx',473;'lambdaEm',520;'na',1.4;'ri',1.518;...
        'sX',s1;'sY',s2;'sZ',s3;...
        'scaleX',res1;'scaleY',res2;'scaleZ',res3});
% ft thereof
kpsf=ft(psf);

%% as an object I want a hollow sphere
% I define it in k-space

kobj=sinc(rr(kpsf)./2).*sinc(rr(kpsf)./0.7);

% this is the sphere
obj=rr(psf).*abs(ift(kobj));
maximum=max(obj);
% in-focus rectangle in right top
obj(83:end,23:43,floor(s3/2))=4*maximum;
% in-focus line on the left
obj(21:21,40:end,floor(s3/2))=12*maximum;

%% shift the object a little bit in z
%kobj=ft(obj);
%dz=0; % shift in pixels -> 1 equals 100nm
%kobj=kobj.*exp(-i*2*pi*zz(kobj,'freq')*dz);
%obj=ift(kobj);

%% excited fluorophores
fluo=obj.*imgrat;

% fluorescence image with structured illumination
strucflimg=ift(ft(fluo).*kpsf);

% widefield fluorescence image
flimg=ift(ft(obj).*kpsf);

% extract focal planes
%iu=real(squeeze(flimg(:,:,floor(s3/2))));
in=real(squeeze(strucflimg(:,:,floor(s3/2))));