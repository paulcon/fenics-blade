%%
% Make the KL bases for the heat flux boundary on the blade.
clear all; close all;

%%
% Make mesh
x = linspace(0,1,31);
z = linspace(0,2,101);
[X,Z] = meshgrid(x,z);
mesh = [X(:) Z(:)];
dlmwrite('xcoord.txt',x,'delimiter',',','precision','%18.16e');
dlmwrite('zcoord.txt',z,'delimiter',',','precision','%18.16e');

%%
% Correlation definition
corr.name='exp';
corr.c0=[0.1 0.01];
corr.c1=0;
corr.sigma=1;

%%
% Get bases
%n = 500;
%[~,KL] = randomfield(corr,mesh,'trunc',n);
%M = bsxfun(@times,KL.bases,KL.sv');
%dlmwrite('klbases.txt',M,'delimiter',',','precision','%18.16e');
%save('klmat.mat','KL');

%%
% 1d bases for testing and dev
corr.c0=0.1;
[~,KL] = randomfield(corr,x');
M = bsxfun(@times,KL.bases,KL.sv');
dlmwrite('klbases_1d.txt',M,'delimiter',',','precision','%18.16e');
save('klmat_1d.mat','KL');