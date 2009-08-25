% simulate clem
%% has to be run as initialization
run './dbl.m'
%% image with grating
grating=mod(rx+3,4)>1;
in=dbl2(grating);
%% widefield image
wf=1;
iu=dbl2(wf);
%% illuminate only in-focus parts of the object
slice_perfect=squeeze(obj(:,:,floor(s3/2))>.1);
iu2=dbl2(slice_perfect);
%% illuminate in-focus parts of the object with a grating
in2=dbl2(slice_perfect & grating);

%%
ft(in)
