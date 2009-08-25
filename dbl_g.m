%% define parameters and global variables here
% use run './dbl_g.m' in the beginning of all scripts that might refer to
% those
res1=100; % nm per pixel
res2=100;
res3=300;

g_phi=3; % pixel, phase shift of grating
s1=128;
s2=128;
s3=128;
rx=xx(s1,s2,'mcorner'); % nm
ry=yy(rx,'mcorner');
Lambda=700; % nm grating period
theta=27*pi/180; % angle of k against x

global grat;
global kgrat;
global kasf0;
global kasf1;
global kasf2;
global psf;
global kpsf;
global obj;
global kobj;