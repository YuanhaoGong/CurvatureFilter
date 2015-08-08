im = imread('lena.png');

Iteration = 10;

tic
[result,energy]=GCFilter(im,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('performance: ', num2str(mytime/Iteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), GCFilter(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Gaussian Curvature Energy'),title('Energy profile')

%% ************************* mean curvature *********************************************
im = imread('lena.png');

Iteration = 10;

tic
[result,energy]=MCFilter(im,Iteration);
mytime = toc;

%% show the running time and the result
mystr = strcat('performance: ', num2str(mytime/Iteration),' seconds per iteration (', num2str(size(im,1)),'X', num2str(size(im,2)), ' image)');
disp(mystr)

figure,imagesc([double(im),result,double(im)-result]), daspect([1,1,1]), colorbar
title('original(left), MCFilter(mid), difference(right)')

figure,plot(energy,'linewidth',4),xlabel('Iteration'), ylabel('Mean Curvature Energy'),title('Energy profile')
