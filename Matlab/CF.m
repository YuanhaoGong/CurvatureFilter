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
im = single(im); im_orig = im; result = im; [m,n]=size(im); Energy = zeros(ItNum,1); 
%% four types of pixels %B = black, W = white, C = circle, T = triangle
[BC_row,BC_col]=meshgrid(2:2:m-1,2:2:n-1);[BT_row,BT_col]=meshgrid(3:2:m-1,3:2:n-1);
[WC_row,WC_col]=meshgrid(2:2:m-1,3:2:n-1);[WT_row,WT_col]=meshgrid(3:2:m-1,2:2:n-1);
BC=sub2ind(size(im),reshape(BC_row,size(BC_row,1)*size(BC_row,2),1), reshape(BC_col,size(BC_col,1)*size(BC_col,2),1));
BT=sub2ind(size(im),reshape(BT_row,size(BT_row,1)*size(BT_row,2),1), reshape(BT_col,size(BT_col,1)*size(BT_col,2),1));
WC=sub2ind(size(im),reshape(WC_row,size(WC_row,1)*size(WC_row,2),1), reshape(WC_col,size(WC_col,1)*size(WC_col,2),1));
WT=sub2ind(size(im),reshape(WT_row,size(WT_row,1)*size(WT_row,2),1), reshape(WT_col,size(WT_col,1)*size(WT_col,2),1));
%% the neighbors' index
BC_pre = BC -1;BC_nex = BC +1;BC_lef = BC -m;BC_rig = BC +m;BC_lu = BC-1-m;BC_ru = BC-1+m;BC_ld=BC+1-m;BC_rd=BC+m+1;
BT_pre = BT -1;BT_nex = BT +1;BT_lef = BT -m;BT_rig = BT +m;BT_lu = BT-1-m;BT_ru = BT-1+m;BT_ld=BT+1-m;BT_rd=BT+m+1;
WC_pre = WC -1;WC_nex = WC +1;WC_lef = WC -m;WC_rig = WC +m;WC_lu = WC-1-m;WC_ru = WC-1+m;WC_ld=WC+1-m;WC_rd=WC+m+1;
WT_pre = WT -1;WT_nex = WT +1;WT_lef = WT -m;WT_rig = WT +m;WT_lu = WT-1-m;WT_ru = WT-1+m;WT_ld=WT+1-m;WT_rd=WT+m+1;
%% Dual Mesh optimization
for i = 1:ItNum
    Energy(i) = mycurv(result); 
    if (i>1) && Energy(i,1) > Energy(i-1,1) % if the energy start to increase
        break;
    end
    result = myfun(result,BC,BC_pre,BC_nex,BC_lef,BC_rig,BC_lu,BC_ld,BC_ru,BC_rd,step);
    result = myfun(result,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd,step);
    result = myfun(result,WC,WC_pre,WC_nex,WC_lef,WC_rig,WC_lu,WC_ld,WC_ru,WC_rd,step);
    result = myfun(result,WT,WT_pre,WT_nex,WT_lef,WT_rig,WT_lu,WT_ld,WT_ru,WT_rd,step);
end
Energy = Energy(1:i,:);
%% %%%%%%%%%%%%%%%%%%%% three projection operaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = proj_TV(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd,step)
res = im; BT5 = 5*im(BT); dist = zeros(size(BT_pre,1),8,'single');
tmp1 = im(BT_pre) + im(BT_nex) - BT5; tmp2 = im(BT_lef) + im(BT_rig) - BT5;
tmp3 = im(BT_pre) + im(BT_lu) + im(BT_ru) - BT5; tmp4 = im(BT_nex) + im(BT_ld) + im(BT_rd) - BT5;
dist(:,1) = tmp1 + im(BT_lu) + im(BT_lef) + im(BT_ld); 
dist(:,2) = tmp1 + im(BT_ru) + im(BT_rig) + im(BT_rd);
dist(:,3) = tmp2 + im(BT_lu) + im(BT_pre) + im(BT_ru); 
dist(:,4) = tmp2 + im(BT_ld) + im(BT_nex) + im(BT_rd); 
dist(:,5) = tmp3 + im(BT_lef) + im(BT_ld); 
dist(:,6) = tmp3 + im(BT_rig) + im(BT_rd);
dist(:,7) = tmp4 + im(BT_lef) + im(BT_lu); 
dist(:,8) = tmp4 + im(BT_rig) + im(BT_ru);
dist = dist/5; %% minimal projection
tmp = abs(dist); [v,ind] = min(tmp,[],2);
tmp = sub2ind(size(dist),(1:size(dist,1))',ind);
tmp = dist(tmp); res(BT) = res(BT) + step*tmp;

function res = proj_MC(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd,step)
res = im; BT8 = 8*im(BT); dist = zeros(size(BT_pre,1),4,'single');
tmp1 = 2.5*(im(BT_pre) + im(BT_nex)) - BT8;
tmp2 = 2.5*(im(BT_lef) + im(BT_rig)) - BT8;
dist(:,1) = tmp1  + 5*im(BT_rig) - im(BT_ru) - im(BT_rd);
dist(:,2) = tmp1  + 5*im(BT_lef) - im(BT_lu) - im(BT_ld);
dist(:,3) = tmp2  + 5*im(BT_pre) - im(BT_lu) - im(BT_ru);
dist(:,4) = tmp2  + 5*im(BT_nex) - im(BT_ld) - im(BT_rd);
dist(:,1:4) = dist(:,1:4)/8; %% minimal projection
tmp = abs(dist); [v,ind] = min(tmp,[],2);
tmp = sub2ind(size(dist),(1:size(dist,1))',ind);
tmp = dist(tmp); res(BT) = res(BT) + step*tmp;

function res = proj_GC(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd,step)
res = im; BT2 = 2*im(BT); BT3 = 3*im(BT);dist = zeros(size(BT_pre,1),8,'single');
dist(:,1) = im(BT_pre) + im(BT_nex) - BT2; dist(:,2) = im(BT_lef) + im(BT_rig) - BT2;
dist(:,3) = im(BT_lu) + im(BT_rd) - BT2; dist(:,4) = im(BT_ld) + im(BT_ru) - BT2;
dist(:,5) = im(BT_pre) + im(BT_lef) + im(BT_lu) - BT3; dist(:,6) = im(BT_pre) + im(BT_rig) + im(BT_ru) - BT3;
dist(:,7) = im(BT_nex) + im(BT_lef) + im(BT_ld) - BT3; dist(:,8) = im(BT_nex) + im(BT_rig) + im(BT_rd) - BT3;
dist(:,1:4) = dist(:,1:4)/2; dist(:,5:8) = dist(:,5:8)/3; %% minimal projection
tmp = abs(dist); [v,ind] = min(tmp,[],2);
tmp = sub2ind(size(dist),(1:size(dist,1))',ind);
tmp = dist(tmp); res(BT) = res(BT) + step*tmp;

function res = proj_BF(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd,step)
res = im; BT2 = 2*im(BT); BT3 = 7*im(BT); dist = zeros(size(BT_pre,1),8,'single');
dist(:,1) = im(BT_pre) + im(BT_nex) - BT2; dist(:,2) = im(BT_lef) + im(BT_rig) - BT2;
dist(:,3) = im(BT_lu) + im(BT_rd) - BT2; dist(:,4) = im(BT_ld) + im(BT_ru) - BT2;
tmp1 = 3*(im(BT_ld) + im(BT_ru)) - BT3; tmp2 = 3*(im(BT_lu) + im(BT_rd)) - BT3;
dist(:,5) = im(BT_pre) + im(BT_lef) - im(BT_lu) + tmp1; 
dist(:,6) = im(BT_pre) + im(BT_rig) - im(BT_ru) + tmp2;
dist(:,7) = im(BT_nex) + im(BT_lef)- im(BT_ld) + tmp2; 
dist(:,8) = im(BT_nex) + im(BT_rig) - im(BT_rd) +tmp1;
dist(:,1:4) = dist(:,1:4)/2; dist(:,5:8) = dist(:,5:8)/7; %% minimal projection
tmp = abs(dist); [v,ind] = min(tmp,[],2);
tmp = sub2ind(size(dist),(1:size(dist,1))',ind);
tmp = dist(tmp); res(BT) = res(BT) + step*tmp;

%% %%%%%%%%%%%%%%%%%%%% three curvature energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function en = curv_TV(im)
t = single(im);[gx,gy]=gradient(t); g = abs(gx) + abs(gy);
en = sum(g(:));
function en = curv_MC(im)
t = single(im);[gx,gy]=gradient(t);[gxx,gxy]=gradient(gx);[gyx,gyy]=gradient(gy);
g = ((1+gy.^2).*gxx + gx.*gy.*(gxy+gyx)+ (1+gx.^2).*gyy)./((1+gx.^2+gy.^2).^1.5)/2;
en = sum(abs(g(:)));
function en = curv_GC(im)
t = single(im);[gx,gy]=gradient(t);[gxx,gxy]=gradient(gx);[gyx,gyy]=gradient(gy);
g = (gxx.*gyy-gxy.*gyx)./((1+gx.^2+gy.^2).^1.5); en = sum(abs(g(:)));

%% %%%%%%%%%%%%%%%%%%%%%%% Fast TV Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result, Energy] = TVFilterFast(im, ItNum, step)
%%% this is the fast TV filter based on box filter to remove the overlapped
%%% computation.
if nargin<3
    step = 1;
end
im = single(im); Energy = zeros(ItNum,1); result = im; [m,n]=size(im);
[col,row]=meshgrid(1:n,1:m); row = reshape(row,m*n,1); col = reshape(col,m*n,1);
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
dist = dist/5; %% minimal projection
tmp = abs(dist); 
[v,ind] = min(tmp,[],3);
ind = reshape(ind,m*n,1);
ind2 = sub2ind(size(dist),row,col,ind);
dm = dist(ind2); 
dm = reshape(dm,m,n);
res = im + step*dm;

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
