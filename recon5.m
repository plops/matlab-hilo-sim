%% start of reconstruction
kin=ft(in);
kiu=ft(iu);
% scale kin and kiu so that ic has no zero order
kappa=sum(abs(kin)./abs(kiu).*(rr(in)<5))./sum(rr(in)<5)
%% correct for the otf
ckin=corr.*kin;
ckiu=corr.*kiu;
%% correlate to find grating positions
ackin=abs(ckin-kappa.*ckiu);
ackiu=abs(ckiu).*(rr(ackiu,'freq')<.16 | abs(xx(ackiu,'freq'))<.06);
%%
kackin=i+newim(512,512);
kackin=kackin-i;
kackiu=kackin;
st=256-64; w=127; en=st+w;
kackin(st:en,st:en)=ft(ackin);
kackiu(st:en,st:en)=ft(ackiu);
kackiu=kackiu;    
%%
kcov=kackin.*conj(kackiu);
coriniu=abs(ift(kcov))
writeim(coriniu*255/max(coriniu),'/home/martin/0807/cov.jpg','JPEG',0);
%%
kcov./abs(kcov) %.*(abs(kcov)>1e-5)
%% suppress zero order by subtracting properly scaled widefield
ic=in-kappa.*iu;
kic=ft(ic)
%% find maximum on the right of the fouriertransform
startx=75;
dic=DampEdge(ic,0.2,2,1,2);
kdic=ft(dic).*corr;
[m,p]=max(abs(kdic(startx:end,:)));
pos=[p(1)+startx,p(2)];
% determine center of mass of the 3x3 region around the maximum
region=abs(kic(pos(1)-1:pos(1)+1,pos(2)-1:pos(2)+1));
region=region-min(region);
cm=[sum(xx(region).*region),sum(yy(region).*region)];
shift=pos+cm-[64,64]
%% multiply by this in object space to shift +1 order into middle
doshift=exp(-i*2*pi*(xx(kic,'freq').*shift(1)+yy(kic,'freq').*shift(2)));
%%
kc=0.045;
r=rr(kic,'freq');
klp=exp(-r.^2/(2*kc^2));
cm=abs(ift(ft(ic.*doshift).*corr.*klp));
ihp=real(ift(ft(iu).*corr.*(1-klp)));
% integrate over ring with radius kc to find eta
ring=abs(ft(besselj(0,2*pi*kc*n.*r)));
ring2=r-1./n<kc & r+1./n>kc;
cring=ring.*ring2;
nring=cring./sum(abs(cring));
eta=sum(abs(ft(ihp))/abs(ft(cm)).*nring);
% combine highpass and lowpass filtered images
ihilo3=2.*eta.*cm+ihp
%% crossection
rothilo=rotation(abs(ft(ihilo3)),40);
sr=size(rothilo);
rothilo(:,sr(2)/2)