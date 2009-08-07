%%
kin=ft(in);
kiu=ft(iu);
kappa=sum(abs(kin)./abs(kiu).*(rr(in)<5))./sum(rr(in)<5)
%%
ic=in-kappa.*iu;
kic=ft(ic)
abs(kic(64,64))
%% find maximum on the right of the fouriertransform
startx=88;
dic=DampEdge(ic,0.2,2,1,2);
kdic=ft(dic);
[m,p]=max(abs(kdic(startx:end,:)));
pos=[p(1)+startx,p(2)];
% determine center of mass of the 3x3 region around the maximum
region=abs(kic(pos(1)-1:pos(1)+1,pos(2)-1:pos(2)+1));
region=region-min(region);
cm=[sum(xx(region).*region),sum(yy(region).*region)];
shift=pos+cm-[64,64]
%%
doshift=exp(-i*2*pi*(xx(kic,'freq').*shift(1)+yy(kic,'freq').*shift(2)));
%%
kc=0.05;
r=rr(kic,'freq');
klp=exp(-r.^2/(2*kc^2));
cm=abs(ift(ft(ic.*doshift).*klp));
ihp=real(ift(ft(iu).*(1-klp)));
ring=abs(ft(besselj(0,2*pi*kc*n.*r)));
ring2=r-1./n<kc & r+1./n>kc;
cring=ring.*ring2;
nring=cring./sum(abs(cring));
eta=sum(abs(ft(ihp))/abs(ft(cm)).*nring);
ihilo3=0.1*eta.*cm+ihp
%% crossection
rothilo=rotation(abs(ft(ihilo3)),40);
sr=size(rothilo);
rothilo(:,sr(2)/2)
%%
ft(ic)
ft(ic-cm)
%%
filter.*kic
%%
filter=gaussf(xx(in,'freq')>0.187,8)
cm=abs(ift(filter.*kic));
isu=cm.*iu