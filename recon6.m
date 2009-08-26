% simulate clem
%% has to be run as initialization
run './dbl.m'
%% image with grating
period=4;
grating=mod(rx,period)>((period/2)-1);
in=dbl2(grating,0);
%% image with shifted grating
gratings=mod(rx+period/2,period)>((period/2)-1);
ins=dbl2(gratings,0);

%% widefield image
wf=1;
iu=dbl2(wf,0);
%% illuminate only in-focus parts of the object
slice_perfect=squeeze(obj(:,:,floor(s3/2))>.1);
iu2=dbl2(slice_perfect,0);
%% illuminate in-focus parts of the object with a grating
in2=dbl2(slice_perfect & grating,0);
%% illuminate in-focus parts of the object with the shifted grating
ins2=dbl2(slice_perfect & gratings,0);

%% project otf along z
skpsf=squeeze(sum(kpsf,[],3));
corr=gaussf((rr(skpsf,'freq')<.42),3)./skpsf; % use this to correct for otf
%%
ft(2.5*in2-iu2).*corr
%% shift right order into the center
%diff=ins2-in2; add=ins2+in2;
%diff=ins-iu; add=iu;
diff=2.5*in2-iu2; add=iu2;
o=ft(ift((ft(diff).*corr)).*exp(-i*2*pi*(xx(in,'freq').*32)));
%%
kc=0.07203;
r=rr(in,'freq');
klp=exp(-r.^2/(2*kc^2)); % low pass filter in k-space
% low resolution slice
ilp=ift(o.*klp)
% high resolution slice
ihp=ift(ft(add).*corr.*(1-klp));

%%
ihilo=ihp+3.63637*ilp
%ihilo(64,:)

%%
1./ft(slice_perfect)