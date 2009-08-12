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
writeimk(log(ackiu).*(ackiu>0),'/home/martin/0807/ackiu.eps');
writeimk(log(ackin).*(ackin>0),'/home/martin/0807/ackin.eps');
%%
kcov=ft(ackin).*conj(ft(ackiu));
kcov_big=newim(512,512)+i-i; % subsample the correlation by zero padding
st=256-64; w=127; en=st+w;
kcov_big(st:en,st:en)=kcov;
cov=abs(ift(kcov_big)); % this contains the cross correlation
%writeim(255*cov./max(cov),'/home/martin/0807/cov.jpg','JPEG',0);
%system('cd /home/martin/0807; convert cov.jpg eps2:cov.eps');
%%
kcov./abs(kcov) %.*(abs(kcov)>1e-5)
%% suppress zero order by subtracting properly scaled widefield
ic=in-kappa.*iu;
kic=ft(ic)
%% find maximum on the right of the correlation
startx=75*4;
[m,p]=max(abs(cov(startx:end,:)));
pos=[p(1)+startx,p(2)];
% determine center of mass of the 3x3 region around the maximum
region=abs(cov(pos(1)-1:pos(1)+1,pos(2)-1:pos(2)+1));
region=region-min(region);
cm=[sum(xx(region).*region),sum(yy(region).*region)];
shift=(pos+cm-[256,256])/4 % divide by 4 to undo subsampling
%% multiply by this in object space to shift +1 order into middle
doshift=exp(-i*2*pi*(xx(ckin,'freq').*shift(1)+yy(ckin,'freq').*shift(2)));
%%
kc=0.052;
r=rr(ckin,'freq');
klp=exp(-r.^2/(2*kc^2)); % low pass filter in k-space
q1=ift(ckin-kappa.*ckiu);
kqs=ft(q1.*doshift);
cm=abs(ift(kqs.*klp));
ihp=real(ift(ft(iu).*corr.*(1-klp)));
% integrate over ring with radius kc to find eta
ring=abs(ft(besselj(0,2*pi*kc*n.*r)));
ring2=r-1./n<kc & r+1./n>kc;
cring=ring.*ring2;
nring=cring./sum(abs(cring));
eta=sum(abs(ft(ihp))/abs(ft(cm)).*nring);
% combine highpass and lowpass filtered images
ihilo3=2.*cm+ihp % don't use eta
%%
writeimk(log(abs(kqs)).*(abs(kqs)>0),'/home/martin/0807/kqs.eps');
writeimk(abs(kqs.*klp),'/home/martin/0807/kqslp.eps');
writeimk(cm,'/home/martin/0807/cm.eps');
writeimk(abs(ft(iu).*corr.*(1-klp)),'/home/martin/0807/kihp3.eps');
writeimk(ihp,'/home/martin/0807/ihp3.eps');
writeimk(ihilo3,'/home/martin/0807/ihilo3.eps');
%% crossection
rothilo=rotation(abs(ft(ihilo3)),40);
sr=size(rothilo);
rothilo(:,sr(2)/2)