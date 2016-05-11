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

im_name = 'lena.png';
%% ************************* Gaussian curvature *********************************************
im = imread(im_name);
if size(im,3)>1
    im = rgb2gray(im);
end

FilterType = 2;
Iteration = 60;

disp('** running time includes the time for computing energy. **')

tic
[result,energy]=CF(im,FilterType, Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('GC filter performance: ', num2str(mytime/size(energy,1)),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure, imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), GCFilter(mid), difference(right)')
figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Gaussian Curvature Energy'),title('Energy profile')

%% ************************* mean curvature *********************************************
im = imread(im_name);
if size(im,3)>1
    im = rgb2gray(im);
end

FilterType = 1;
Iteration = 60;

tic
[result,energy]=CF(im,FilterType,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('MC filter performance: ', num2str(mytime/size(energy,1)),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure, imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), MCFilter(mid), difference(right)')
figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Mean Curvature Energy'),title('Energy profile')

%% ************************* Bernstein Filter also minimizes mean curvature *********************************************
im = imread(im_name);
if size(im,3)>1
    im = rgb2gray(im);
end

FilterType = 4;
Iteration = 60;

tic
[result,energy]=CF(im, FilterType, Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('BF filter performance: ', num2str(mytime/size(energy,1)),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure, imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), BFilter(mid), difference(right)')
figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Mean Curvature Energy'),title('Energy profile')

%% ************************* TV Filter minimizes Total Variation *********************************************
im = imread(im_name);
if size(im,3)>1
    im = rgb2gray(im);
end

FilterType = 0;
Iteration = 60;

tic
[result,energy]=CF(im, FilterType, Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('TV filter performance: ', num2str(mytime/size(energy,1)),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure, imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), TVFilter(mid), difference(right)')
figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('TV Energy'),title('Energy profile')


%% ************************************************************************
%
%              Filter Solver for Variational Models
%
% (these filters only minimize the regularization term. When they are used as
% solver for generic variational models, data fitting term has to be considered.)
%
%             one instance of Algorithm 17 in my PhD thesis
%% ************************************************************************

im = imread(im_name);
if size(im,3)>1
    im = rgb2gray(im);
end

MaxIteration = 60;
DataFitOrder = 1.3; %fractional order
Lambda = 3;
FilterType = 0;

tic
[result,energy]=Solver(im, FilterType, DataFitOrder, Lambda, MaxIteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('solver performance: ', num2str(mytime/size(energy,1)),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), Result(mid), difference(right)')

x=0:size(energy,1)-1; x=x';
figure,plot(x,energy(:,1),'linewidth',4),xlabel('Iteration'), ylabel('Energy'),title('Energy profile'),hold on
plot(x,energy(:,2),'linewidth',4),plot(x,energy(:,3),'linewidth',4), legend('Total Energy','DataFit Energy','Regularization Energy','location','west')
legend('boxoff')
