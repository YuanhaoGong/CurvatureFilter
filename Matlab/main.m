% =========================================================================
% 
%                           Curvature Filter 
% 
% *************************************************************************  
% 
%         @phdthesis{gong:phd, 
%          title={Spectrally regularized surfaces}, 
%          author={Gong, Yuanhao}, 
%          year={2015}, 
%          school={ETH Zurich, Nr. 22616},
%          note={http://dx.doi.org/10.3929/ethz-a-010438292}}
% 
% =========================================================================

% this demo shows four edge-preserving filters and how to use them solve variational models 

%% ************************* Gaussian curvature *********************************************
im = imread('lena.png');

Iteration = 60;

disp('** running time includes the time for computing energy. **')

tic
[result,energy]=GCFilter(im,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('GC filter performance: ', num2str(mytime/Iteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), GCFilter(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Gaussian Curvature Energy'),title('Energy profile')

%% ************************* mean curvature *********************************************
im = imread('lena.png');

Iteration = 60;

tic
[result,energy]=MCFilter(im,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('MC filter performance: ', num2str(mytime/Iteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), MCFilter(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Mean Curvature Energy'),title('Energy profile')

%% ************************* Bernstein Filter also minimizes mean curvature *********************************************
im = imread('lena.png');

Iteration = 60;

tic
[result,energy]=BernsteinFilter(im,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('BF filter performance: ', num2str(mytime/Iteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), BFilter(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Mean Curvature Energy'),title('Energy profile')

%% ************************* TV Filter minimizes Total Variation *********************************************
im = imread('lena.png');

Iteration = 60;

tic
[result,energy]=TVFilter(im,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('TV filter performance: ', num2str(mytime/Iteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), TVFilter(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('TV Energy'),title('Energy profile')




%% ************************************************************************
%
%              Local Filter Solver for Variational Models
%
%% ************************************************************************

im = imread('lena.png');

MaxIteration = 60;
DataFitOrder = 1.3; %fractional order
Lambda = 3;
FilterType = 0;

tic
[result,energy]=Solver(im, FilterType, DataFitOrder, Lambda, MaxIteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('solver performance: ', num2str(mytime/MaxIteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), Result(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Total Energy'),title('Energy profile')
