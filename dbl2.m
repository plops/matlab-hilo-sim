% when dbl.m has been run once use this function to generate images with
% the shifted grating. this function relies on results that already have
% been computed by dbl.m
function in=dbl2(g_phi_grad)
    run './dbl_g.m';
    g=1+.5*cos(2*pi/Lambda*(cos(theta)*rx*res1+sin(theta)*ry*res2)+g_phi_grad*pi/180);
    
    % update illumination
    grat(:,:,floor(s3/2))=g(:,:);
    kgrat=ft(grat);
    imgratx=ift(kgrat.*kasf0);
    imgraty=ift(kgrat.*kasf1);
    imgratz=ift(kgrat.*kasf2);
    imgrat=abs(imgratx).^2+abs(imgraty).^2+abs(imgratz).^2;
    %% excited fluorophores
    fluo=obj.*imgrat;
    % the structured illuminated stuff in image space
    strucflimg=ift(ft(fluo).*kpsf);
    % the result
    in=real(squeeze(strucflimg(:,:,floor(s3/2))));
end