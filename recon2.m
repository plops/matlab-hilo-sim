%%
si=size(in);
n=si(1);
r=rr(in,'freq');
ratio=in./iu;
ft(ratio)

%%
filter=gaussf(xx(in,'freq')>0.2,8);
cm=abs(ift(filter.*ft(ratio)));
isu=cm.*iu

%%
kc=1/(2*pi);
klp=exp(-r.^2/(2*kc^2));
ilp=real(ift(klp.*ft(isu)))
%%
ihp=real(ift((1-klp).*ft(iu)))
%%
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
ihilo=0.7*eta.*ilp+ihp
%%
ft(ihilo)