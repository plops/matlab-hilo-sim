%%
si=size(in);
n=si(1);
r=rr(in,'freq');
ratio=in./iu;
ft(ratio)

%% select R+ and construct isu
filter=gaussf((xx(in,'freq')>0.1) & (xx(in,'freq')<0.27),4);
ftratio=ft(ratio);
cm=abs(ift(filter.*ftratio));
isu=cm.*iu

%% calculate low pass filtered low-res section
kc=.07;
klp=exp(-r.^2/(2*kc^2));
ilp=real(ift(klp.*ft(isu)))
%% calculate high pass filtered high-res section
ihp=real(ift((1-klp).*ft(iu)))
%% construct the circle for integration in k-space
ring=abs(ft(besselj(0,2*pi*kc*n.*r)));
ring2=r-1./n<kc & r+1./n>kc;
cring=ring.*ring2;
nring=cring./sum(abs(cring))
%%
%eta=sum(abs(ft(ihp)).*nring)./sum(abs(ft(ilp)).*nring)
%%
%eta=sum(abs(ft(ihp)).*nring./((abs(ft(ilp)).*nring)+(nring==0)))
%% i think that is the right way
eta=sum(abs(ft(ihp))/abs(ft(ilp)).*nring)
%%
ihilo=eta.*ilp+ihp
%%
ft(ihilo)
%% scalebar
sbx=100;
sby=110;
sbh=3;
sbw=20; % 2mu
objsb=real(ilp);
%log(abs(filter.*ftratio)+.000001);
% real(ratio)
objsb=(objsb-min(objsb))*255/(max(objsb)-min(objsb));
objsb(sbx:sbx+sbw,sby:sby+sbh)=255;
writeim(objsb,'/home/martin/0807/ilp2.eps','EPS',0);