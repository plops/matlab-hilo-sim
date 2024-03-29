%%
si=size(in);
n=si(1);
kc=0.04;
r=rr(in,'freq');
klp=exp(-r.^2/(2*kc^2));
khp=1-klp;
ina=ift(ft(in).*klp);
in2a=ift(ft(in.^2).*klp);
iua=ift(ft(iu).*klp);
iu2a=ift(ft(iu.^2).*klp);
cs2=in2a.*iua.^2./(ina.^2.*iu2a)-1;
isu=sqrt(cs2).*iua;
kihp=ft(iu).*khp;
kilp=ft(isu).*klp;
ring=real(ft(besselj(0,2*pi*kc*n.*r)));
ring2=r-1./n<kc & r+1./n>kc;
cring=ring.*ring2;
nring=cring./sum(cring); % normalized ring with radius kc
%hring=abs(cring.*kihp);
%lring=abs(cring.*kilp);
%eta=sum(hring./lring)./sum(cring);
%lr=sum(abs(cring.*kilp));
%ur=sum(abs(cring.*kihp));
%eta=ur./lr;
%eta=1;
eta=sum(abs(kihp)./abs(kilp).*nring)
kihilo=eta.*kilp+kihp;
ihilo=ift(kihilo)
mix(:,64);
%% scalebar
sbx=100;
sby=110;
sbh=3;
sbw=20; % 2mu
objsb=real(ihilo); %real(ift(kihp));
%log(abs(ihilo+nring/100)); %sqrt(abs(cs2)) abs(isu) (ift(kilp)) abs(ift(kilp))
objsb=(objsb-min(objsb))*255/(max(objsb)-min(objsb));
objsb(sbx:sbx+sbw,sby:sby+sbh)=255;
writeim(objsb,'/home/martin/0807/ihilo.eps','EPS',0);
%%
