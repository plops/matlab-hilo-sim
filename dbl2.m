% when dbl.m has been run once use this function to generate images with
% different illumination this function relies on results that already have
% been computed by dbl.m
function in=dbl2(illum)
    run './dbl_g.m';
    % update illumination
    grat(:,:,floor(s3/2))=illum(:,:);
    kgrat=ft(grat);
    imgratx=ift(kgrat.*kasf0);
    imgraty=ift(kgrat.*kasf1);
    imgratz=ift(kgrat.*kasf2);
    imgrat=abs(imgratx).^2+abs(imgraty).^2+abs(imgratz).^2;
    %% excited fluorophores
    fluo=obj.*imgrat;
    %fluo=obj.*grat;
    % the structured illuminated stuff in image space
    strucflimg=ift(ft(fluo).*kpsf);
    % the result
    in=real(squeeze(strucflimg(:,:,floor(s3/2))));
end