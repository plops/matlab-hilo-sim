% simulate clem
%% has to be run as initialization
run './dbl.m'
laser=1;
%% image with grating
period=4;
grating=mod(rx,period)>((period/2)-1);
in=dbl2(grating,laser);
%% image with shifted grating
gratings=mod(rx+period/2,period)>((period/2)-1);
ins=dbl2(gratings,laser);

%% widefield image
wf=1;
iu=dbl2(wf,laser);
%% illuminate only in-focus parts of the object
%gaussf
slice_perfect=bdilation(squeeze(obj(:,:,floor(s3/2))>.1),1);
iu2=dbl2(slice_perfect,laser);
%% illuminate in-focus parts of the object with a grating
in2=dbl2(slice_perfect & grating,laser);
%% illuminate in-focus parts of the object with the shifted grating
ins2=dbl2(slice_perfect & gratings,laser);

%% for the clem image remove in-focus rectangle in right top
rectangle=newim(slice_perfect)>1;
rectangle(83:end,23:43)=1;
slice_clem=slice_perfect & ~ bdilation(rectangle);
iu3=dbl2(slice_clem,laser);
%%
in3=dbl2(slice_clem & grating,laser);

%% only illuminate the bright rectangle
not_slice_clem=bdilation(rectangle);
iu4=dbl2(not_slice_clem,laser);
%%
in4=dbl2(not_slice_clem & grating,laser);

%% project otf along z
skpsf=squeeze(sum(kpsf,[],3));
correct=gaussf((rr(skpsf,'freq')<.42),3)./skpsf; % use this to correct for otf
%%
ft(2.06612*in2-iu2).*corr
%%
ft(2*in3-iu3).*corr
%% only in-focus
ihilo2=hilo_combine(2*in2-iu2,iu2)
%% leave out rectangle
ihilo3=hilo_combine(2*in3-iu3,iu3)
%% only rectangle
ihilo4=hilo_combine(2*in4-iu4,iu4)
%%
ihilo3+ihilo4 %-ihilo2
%% shift right order into the center
%diff=ins2-in2; add=ins2+in2;
%diff=ins-iu; add=iu;
%diff=2.06612*in2-iu2; add=iu2;
diff=2.*in3-iu3; add=iu3;
%diff=2*in-iu; add=iu;
o=ft(ift((ft(diff).*corr)).*exp(-i*2*pi*(xx(in,'freq').*32)));
%%
kc=0.06;
r=rr(in,'freq');
klp=exp(-r.^2/(2*kc^2)); % low pass filter in k-space
% low resolution slice
ilp=abs(ift(o.*klp))
% high resolution slice
ihp=abs(ift(ft(add).*corr.*(1-klp)));

%%
ihilo=ihp+6.90454*ilp
ihilo(64,:)

%%
1./ft(slice_perfect)