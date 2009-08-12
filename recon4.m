% I thought it would be possible to divide by the highpass filtered 
% grating in order to move the modulated part of the image into DC.
% However, this doesn't work as good as I hoped. There are considerable
% residual peaks left.

%% start of reconstruction
kin=ft(in);
kiu=ft(iu);
% scale kin and kiu
kappa=sum(abs(kin)./abs(kiu).*(rr(in)<5))./sum(rr(in)<5)

%% correct for the otf
ckin=corr.*kin;
ckiu=corr.*kiu;

%% remove zero order from in
q=in-kappa.*iu;
kq=ft(q);

%%
r=rr(kic,'freq');
kc=.07;
klp=exp(-r.^2/(2*kc^2))

%% extract illumination in focal plane
i0=squeeze(imgrat(:,:,128));
ii=ift((1-klp).*ft(i0))

%%
qq=ift(kq.*corr);
o0=qq./ii;
%%
o1=ift(ft(o0).*exp(-r.^2/(2*0.05^2)))