function ihilo=hilo_combine(diff,add)
    run './dbl_g.m';
    %diff=2*in-iu; add=iu;  
    o=ft(ift((ft(diff).*correct)).*exp(-i*2*pi*(xx(diff,'freq').*32)));
    %
    kc=0.06;
    r=rr(diff,'freq');
    klp=exp(-r.^2/(2*kc^2)); % low pass filter in k-space
    % low resolution slice
    ilp=abs(ift(o.*klp));
    % high resolution slice
    ihp=abs(ift(ft(add).*correct.*(1-klp)));

    %
    ihilo=ihp+6.90454*ilp;
end