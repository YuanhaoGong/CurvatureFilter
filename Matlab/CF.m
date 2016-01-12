function [result, Energy] = CF(im, FilterType, ItNum, step)
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
% FilterType: 0(Total Variation), 1(Mean Curvature), 2(Gaussian Curvature)
%             4(Bernstein Filter), 5(Fast TV Filter)
if (nargin~=3) && (nargin~=4)
    disp('Input are not correct.'), return;
end
if nargin==3
    step = 1;
end
if (size(im,3)>1)
    disp('This only works for gray scale image.'), return;
end
switch FilterType
    case 0
        myfun = @proj_TV; mycurv = @curv_TV;
    case 1
        myfun = @proj_MC; mycurv = @curv_MC;
    case 2
        myfun = @proj_GC; mycurv = @curv_GC;
    case 4
        myfun = @proj_BF; mycurv = @curv_MC;
    case 5
        [result, Energy] = TVFilterFast(im, ItNum, step); return;
   otherwise
      disp('Filter Type is not correct.'), return;
end
%pad to even size
orig = im; [orig_r, orig_c]=size(orig); m=ceil(orig_r/2)*2; n=ceil(orig_c/2)*2;im = zeros(m,n,'single'); 
im(1:orig_r,1:orig_c)=single(orig); im(m,:) = im(m-1,:); im(:,n) = im(:,n-1); result = im; Energy = zeros(ItNum,1); 

%% four types of pixels %B = black, W = white, C = circle, T = triangle; r and c indicate row and column
BC_r = uint32(2:2:m-2); BC_c = uint32(2:2:n-2); BT_r = uint32(3:2:m-1); BT_c = uint32(3:2:n-1);
WC_r = uint32(2:2:m-2); WC_c = uint32(3:2:n-1); WT_r = uint32(3:2:m-1); WT_c = uint32(2:2:n-2);
%% the neighbors' index
BC_pre = BC_r-1; BC_nex = BC_r+1; BC_lef = BC_c-1; BC_rig = BC_c+1;
BT_pre = BT_r-1; BT_nex = BT_r+1; BT_lef = BT_c-1; BT_rig = BT_c+1;
WC_pre = WC_r-1; WC_nex = WC_r+1; WC_lef = WC_c-1; WC_rig = WC_c+1;
WT_pre = WT_r-1; WT_nex = WT_r+1; WT_lef = WT_c-1; WT_rig = WT_c+1;
[row,col] = ndgrid(1:size(BT_r,2),1:size(BT_c,2)); row=uint32(row); col=uint32(col);
%% Dual Mesh optimization
for i = 1:ItNum
    Energy(i) = mycurv(result); 
    if (i>1) && Energy(i,1) > Energy(i-1,1) % if the energy start to increase
        break;
    end
    result = myfun(result,BC_r,BC_c,BC_pre,BC_nex,BC_lef,BC_rig,row,col,step);
    result = myfun(result,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step);
    result = myfun(result,WC_r,WC_c,WC_pre,WC_nex,WC_lef,WC_rig,row,col,step);
    result = myfun(result,WT_r,WT_c,WT_pre,WT_nex,WT_lef,WT_rig,row,col,step);
end
%unpad
Energy = Energy(1:i,:); result = result(1:orig_r,1:orig_c);
%% %%%%%%%%%%%%%%%%%%%% three projection operaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = proj_TV(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT5 = 5*im(BT_r,BT_c); dist = zeros([size(BT5),8],'single');
tmp1 = im(BT_pre,BT_c) + im(BT_nex,BT_c) - BT5; tmp2 = im(BT_r,BT_lef) + im(BT_r,BT_rig) - BT5;
tmp3 = im(BT_pre,BT_c) + im(BT_pre,BT_lef) + im(BT_pre,BT_rig) - BT5; 
tmp4 = im(BT_nex,BT_c) + im(BT_nex,BT_lef) + im(BT_nex,BT_rig) - BT5;

dist(:,:,1) = tmp1 + im(BT_pre,BT_lef) + im(BT_r,BT_lef) + im(BT_nex,BT_lef); 
dist(:,:,2) = tmp1 + im(BT_pre,BT_rig) + im(BT_r,BT_rig) + im(BT_nex,BT_rig);
dist(:,:,3) = tmp2 + im(BT_pre,BT_lef) + im(BT_pre,BT_c) + im(BT_pre,BT_rig); 
dist(:,:,4) = tmp2 + im(BT_nex,BT_lef) + im(BT_nex,BT_c) + im(BT_nex,BT_rig); 
dist(:,:,5) = tmp3 + im(BT_r,BT_lef) + im(BT_nex,BT_lef); 
dist(:,:,6) = tmp3 + im(BT_r,BT_rig) + im(BT_nex,BT_rig);
dist(:,:,7) = tmp4 + im(BT_r,BT_lef) + im(BT_pre,BT_lef); 
dist(:,:,8) = tmp4 + im(BT_r,BT_rig) + im(BT_pre,BT_rig);
%% minimal projection
tmp = abs(dist); [v,ind] = min(tmp,[],3);
index = sub2ind(size(dist), row, col, ind);
dm = step/5*dist(index); 
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

function res = proj_MC(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT8 = 8*im(BT_r,BT_c); dist = zeros([size(BT8),4],'single');
tmp1 = 2.5*(im(BT_pre,BT_c) + im(BT_nex,BT_c)) - BT8;
tmp2 = 2.5*(im(BT_r,BT_lef) + im(BT_r,BT_rig)) - BT8;
dist(:,:,1) = tmp1  + 5*im(BT_r,BT_rig) - im(BT_pre,BT_rig) - im(BT_nex,BT_rig);
dist(:,:,2) = tmp1  + 5*im(BT_r,BT_lef) - im(BT_pre,BT_lef) - im(BT_nex,BT_lef);
dist(:,:,3) = tmp2  + 5*im(BT_pre,BT_c) - im(BT_pre,BT_lef) - im(BT_pre,BT_rig);
dist(:,:,4) = tmp2  + 5*im(BT_nex,BT_c) - im(BT_nex,BT_lef) - im(BT_nex,BT_rig);
tmp = abs(dist); [v,ind] = min(tmp,[],3);
index = sub2ind(size(dist),row,col,ind);
dm = step/8*dist(index); 
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

function res = proj_GC(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT2 = 2*im(BT_r,BT_c); BT3 = 3*im(BT_r,BT_c);dist = zeros([size(BT2),8],'single');
dist(:,:,1) = im(BT_pre,BT_c) + im(BT_nex,BT_c) - BT2; 
dist(:,:,2) = im(BT_r,BT_lef) + im(BT_r,BT_rig) - BT2;
dist(:,:,3) = im(BT_pre,BT_lef) + im(BT_nex,BT_rig) - BT2; 
dist(:,:,4) = im(BT_nex,BT_lef) + im(BT_pre,BT_rig) - BT2;
dist(:,:,5) = im(BT_pre,BT_c) + im(BT_r,BT_lef) + im(BT_pre,BT_lef) - BT3; 
dist(:,:,6) = im(BT_pre,BT_c) + im(BT_r,BT_rig) + im(BT_pre,BT_rig) - BT3;
dist(:,:,7) = im(BT_nex,BT_c) + im(BT_r,BT_lef) + im(BT_nex,BT_lef) - BT3; 
dist(:,:,8) = im(BT_nex,BT_c) + im(BT_r,BT_rig) + im(BT_nex,BT_rig) - BT3;
dist(:,:,1:4) = dist(:,:,1:4)*1.5; %% scale to the same level
tmp = abs(dist); 
[v,ind] = min(tmp,[],3); 
index = sub2ind(size(dist),row,col,ind);
dm = single(step/3)*dist(index); 
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

function res = proj_BF(im,BT_r,BT_c,BT_pre,BT_nex,BT_lef,BT_rig,row,col,step)
res = im; BT2 = 2*im(BT_r,BT_c); BT7 = 7*im(BT_r,BT_c); 
dist = zeros([size(BT2),6],'single');
dist(:,:,1) = im(BT_pre,BT_c) + im(BT_nex,BT_c) - BT2; 
dist(:,:,2) = im(BT_r,BT_lef) + im(BT_r,BT_rig) - BT2;
tmp1 = 3*(im(BT_nex,BT_lef) + im(BT_pre,BT_rig)) - BT7;
tmp2 = 3*(im(BT_pre,BT_lef) + im(BT_nex,BT_rig)) - BT7;
dist(:,:,3) = im(BT_pre,BT_c) + im(BT_r,BT_lef) - im(BT_pre,BT_lef) + tmp1; 
dist(:,:,4) = im(BT_pre,BT_c) + im(BT_r,BT_rig) - im(BT_pre,BT_rig) + tmp2;
dist(:,:,5) = im(BT_nex,BT_c) + im(BT_r,BT_lef)- im(BT_nex,BT_lef) + tmp2; 
dist(:,:,6) = im(BT_nex,BT_c) + im(BT_r,BT_rig) - im(BT_nex,BT_rig) +tmp1;
dist(:,:,1:2) = 10/3*dist(:,:,1:2); %% scale to the same level
tmp = abs(dist); [v,ind] = min(tmp,[],3);
index = sub2ind(size(dist),row,col,ind);
dm = step/10*dist(index); 
res(BT_r,BT_c) = res(BT_r,BT_c) + dm;

%% %%%%%%%%%%%%%%%%%%%% three curvature energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function en = curv_TV(im)
[gx,gy]=gradient(im); g = abs(gx) + abs(gy);
en = sum(g(:));
function en = curv_MC(im)
[gx,gy]=gradient(im);[gxx,gxy]=gradient(gx);[gyx,gyy]=gradient(gy);
g = ((1+gy.^2).*gxx + gx.*gy.*(gxy+gyx)+ (1+gx.^2).*gyy)./((1+gx.^2+gy.^2).^1.5)/2;
en = sum(abs(g(:)));
function en = curv_GC(im)
[gx,gy]=gradient(im);[gxx,gxy]=gradient(gx);[gyx,gyy]=gradient(gy);
g = (gxx.*gyy-gxy.*gyx)./((1+gx.^2+gy.^2).^1.5); en = sum(abs(g(:)));

%% %%%%%%%%%%%%%%%%%%%%%%% Fast TV Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, Energy] = TVFilterFast(im, ItNum, step)
%%% this is the fast TV filter based on box filter to remove the overlapped
%%% computation.
if nargin<3
    step = 1;
end
im = single(im); Energy = zeros(ItNum,1); result = im; [m,n]=size(im);
[row,col]=ndgrid(1:m,1:n);
%% not dual Mesh optimization
for it = 1:ItNum
    Energy(it) = curv_TV(result);
    [total, vert, horiz, diag] = Half_Box(result);
    result = SimpleUpdate(result, total, vert, horiz, diag, step, row, col);
end
function res = SimpleUpdate(im, total, vert, horiz, diag, step, row, col)
res = im; BT5 = total - 5*im; BT6 = total - 6*im; [m,n]=size(im); dist = zeros(m,n,8,'single');
dist(1:m-1,:,1) = BT6(1:m-1,:) - horiz(2:m,:); 
dist(2:m,:,2) = BT6(2:m,:) - horiz(1:m-1,:); 
dist(:,2:n,3) = BT6(:,2:n) - vert(:,1:n-1); 
dist(:,1:n-1,4) = BT6(:,1:n-1) - vert(:,2:n); 
dist(:,:,5) = BT5 - diag; 
dist(:,1:n-1,6) = BT5(:,1:n-1) - diag(:,2:n); 
dist(1:m-1,1:n-1,7) = BT5(1:m-1,1:n-1) - diag(2:m,2:n) ; 
dist(1:m-1,:,8) = BT5(1:m-1,:) - diag(2:m,:); 
%% minimal projection
tmp = abs(dist); 
[v,ind] = min(tmp,[],3);
ind2 = sub2ind(size(dist),row,col,ind);
dm = step/5*dist(ind2); 
res = im + dm;

function [total, vert, horiz, diag] = Half_Box(im)
%compute the total, vertical and horizontal sum in a 3X3 window
%compute the diag(left up corner) sum
[m,n]=size(im);vert = zeros(m,n,'single'); horiz = vert; diag = vert;
total = myboxfilter(im);
imCum = cumsum(im,1);
vert(3:m-1,:) = imCum(4:m,:) - imCum(1:m-3,:);
imCum = cumsum(im,2);
horiz(:,3:n-1) = imCum(:,4:n) - imCum(:,1:n-3);
diag(2:m,2:n) = im(2:m,2:n) + im(1:m-1,2:n) + im(2:m,1:n-1) + im(1:m-1,1:n-1);

function imDst = myboxfilter(imSrc)
[hei, wid] = size(imSrc); imDst = zeros(size(imSrc),'single'); r=1;
%cumulative sum over Y axis
imCum = cumsum(imSrc, 1);
%difference over Y axis
imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);
%cumulative sum over X axis
imCum = cumsum(imDst, 2);
%difference over Y axis
imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
