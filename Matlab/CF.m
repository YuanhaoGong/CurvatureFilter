function [result, Energy] = CF(im, FilterType, ItNum, stepsize)
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
%             4(Bernstein Filter), 5(Half Laplace Filter)
if (nargout>1)
    if_record_energy = true;
else
    if_record_energy = false;
end
switch nargin
    case 2
        ItNum = 10; stepsize = 1;
    case 3
        stepsize = 1;
    case 4
    otherwise
        disp('Input parameters are not correct.'), result=0; Energy=0; return;
end

switch FilterType
    case 0
        myfun = @proj_TV; mycurv = @curv_TV;
    case 1
        myfun = @proj_MC; mycurv = @curv_MC;
    case 2
        myfun = @proj_GC; mycurv = @curv_GC;
    case 3
        myfun = @proj_HL; mycurv = @curv_MC;
   otherwise
      disp('Filter Type is not correct.'), return;
end
%pad to even size
orig = single(im); [orig_r, orig_c, orig_z]=size(orig); m=ceil(orig_r/2)*2; n=ceil(orig_c/2)*2;
im = zeros(m,n,orig_z,'single'); im(1:orig_r,1:orig_c,1:orig_z)=orig; 
im(m,:,:) = im(m-1,:,:); im(:,n,:) = im(:,n-1,:); 
%init
result = im; Energy = zeros(ItNum,orig_z);
%% row and column subscript for four types of pixels 
%B = black, W = white, C = circle, T = triangle; r and c indicate row and column
BC_r = uint32(2:2:m-2); BC_c = uint32(2:2:n-2); BT_r = uint32(3:2:m-1); BT_c = uint32(3:2:n-1);
WC_r = uint32(2:2:m-2); WC_c = uint32(3:2:n-1); WT_r = uint32(3:2:m-1); WT_c = uint32(2:2:n-2);

[row,col] = ndgrid(1:size(BT_r,2),1:size(BT_c,2)); row=uint32(row); col=uint32(col);
%% Dual Mesh optimization
for ch = 1:orig_z
    for i = 1:ItNum
        if(if_record_energy) Energy(i, ch) = mycurv(result); %record curvature energy
            if (i>1) && Energy(i,ch) > Energy(i-1,ch) % if the energy start to increase
                break;
            end
        end
        result(:,:,ch) = myfun(result(:,:,ch),BC_r,BC_c,row,col,stepsize);
        result(:,:,ch) = myfun(result(:,:,ch),BT_r,BT_c,row,col,stepsize);
        result(:,:,ch) = myfun(result(:,:,ch),WC_r,WC_c,row,col,stepsize);
        result(:,:,ch) = myfun(result(:,:,ch),WT_r,WT_c,row,col,stepsize);
    end
end
%unpad
result = result(1:orig_r,1:orig_c,1:orig_z);
%% %%%%%%%%%%%%%%%%%%%% three projection operaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = proj_TV(im,BT_r,BT_c,row,col,step)
res = im; 
BT = im(BT_r,BT_c);
BT5 = 5*BT; 
dist = zeros([size(BT),8],'single');
%eight neighbors
im_pre_c = im(BT_r-1,BT_c);
im_nex_c = im(BT_r+1,BT_c);
im_r_lef = im(BT_r,BT_c-1);
im_r_rig = im(BT_r,BT_c+1);
im_pre_lef = im(BT_r-1,BT_c-1);
im_pre_rig = im(BT_r-1,BT_c+1);
im_nex_lef = im(BT_r+1,BT_c-1);
im_nex_rig = im(BT_r+1,BT_c+1);
%common
tmp1 = im_pre_c + im_nex_c - BT5; 
tmp2 = im_r_lef + im_r_rig - BT5;
tmp3 = im_pre_c + im_pre_lef + im_pre_rig - BT5; 
tmp4 = im_nex_c + im_nex_lef + im_nex_rig - BT5;
%compute all possible projection distances
dist(:,:,1) = tmp1 + im_pre_lef + im_r_lef + im_nex_lef;
dist(:,:,2) = tmp1 + im_pre_rig + im_r_rig + im_nex_rig;
dist(:,:,3) = tmp2 + im_pre_lef + im_pre_c + im_pre_rig; 
dist(:,:,4) = tmp2 + im_nex_lef + im_nex_c + im_nex_rig; 
dist(:,:,5) = tmp3 + im_r_lef + im_nex_lef; 
dist(:,:,6) = tmp3 + im_r_rig + im_nex_rig;
dist(:,:,7) = tmp4 + im_r_lef + im_pre_lef; 
dist(:,:,8) = tmp4 + im_r_rig + im_pre_rig;
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[~,ind] = min(tmp,[],3);
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + dim1.*uint32(col-1) + dim2.*uint32(ind-1);
dm = single(step/5)*dist(index); 
%update current pixels
res(BT_r,BT_c) = BT + dm;

function res = proj_MC(im,BT_r,BT_c,row,col,step)
res = im; 
BT = im(BT_r,BT_c);
BT8 = 8*BT; 
dist = zeros([size(BT),4],'single');
%eight neighbors
im_pre_c = im(BT_r-1,BT_c);
im_nex_c = im(BT_r+1,BT_c);
im_r_lef = im(BT_r,BT_c-1);
im_r_rig = im(BT_r,BT_c+1);
im_pre_lef = im(BT_r-1,BT_c-1);
im_pre_rig = im(BT_r-1,BT_c+1);
im_nex_lef = im(BT_r+1,BT_c-1);
im_nex_rig = im(BT_r+1,BT_c+1);
%common
tmp1 = 2.5*(im_pre_c + im_nex_c) - BT8;
tmp2 = 2.5*(im_r_lef + im_r_rig) - BT8;
%compute all possible projection distances
dist(:,:,1) = tmp1  + 5*im_r_rig - im_pre_rig - im_nex_rig;
dist(:,:,2) = tmp1  + 5*im_r_lef - im_pre_lef - im_nex_lef;
dist(:,:,3) = tmp2  + 5*im_pre_c - im_pre_lef - im_pre_rig;
dist(:,:,4) = tmp2  + 5*im_nex_c - im_nex_lef - im_nex_rig;
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[~,ind] = min(tmp,[],3);
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + dim1.*uint32(col-1) + dim2.*uint32(ind-1);
dm = single(step/8)*dist(index); 
%update current pixels
res(BT_r,BT_c) = BT + dm;

function res = proj_HL(im,BT_r,BT_c,row,col,step)
res = im; 
BT = im(BT_r,BT_c);
BT6 = 6*BT; 
dist = zeros([size(BT),4],'single');
%four neighbors
im_leftup = im(BT_r-1,BT_c-1)+im(BT_r-1,BT_c)+im(BT_r,BT_c-1);
im_leftdn = im(BT_r+1,BT_c-1)+im(BT_r+1,BT_c)+im(BT_r,BT_c-1);
im_rightup = im(BT_r-1,BT_c+1)+im(BT_r-1,BT_c)+im(BT_r,BT_c+1);
im_rightdn = im(BT_r+1,BT_c+1)+im(BT_r+1,BT_c)+im(BT_r,BT_c+1);
%compute all possible projection distances
dist(:,:,1) = im_leftup  + im_leftdn - BT6;
dist(:,:,2) = im_rightup  + im_rightdn - BT6;
dist(:,:,3) = im_leftup  + im_rightup - BT6;
dist(:,:,4) = im_leftdn  + im_rightdn - BT6;
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[~,ind] = min(tmp,[],3);
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + dim1.*uint32(col-1) + dim2.*uint32(ind-1);
dm = single(step/6)*dist(index); 
%update current pixels
res(BT_r,BT_c) = BT + dm;

function res = proj_GC(im,BT_r,BT_c,row,col,step)
res = im; 
BT = im(BT_r,BT_c);
BT2 = 2*BT; 
BT3 = 1.5*BT2; 
dist = zeros([size(BT),8],'single');
%eight neighbors
im_pre_c = im(BT_r-1,BT_c);
im_nex_c = im(BT_r+1,BT_c);
im_r_lef = im(BT_r,BT_c-1);
im_r_rig = im(BT_r,BT_c+1);
im_pre_lef = im(BT_r-1,BT_c-1);
im_pre_rig = im(BT_r-1,BT_c+1);
im_nex_lef = im(BT_r+1,BT_c-1);
im_nex_rig = im(BT_r+1,BT_c+1);
%common
tmp1 = im_pre_c - BT3;
tmp2 = im_nex_c - BT3;
dist(:,:,1) = im_pre_c + im_nex_c - BT2; 
dist(:,:,2) = im_r_lef + im_r_rig - BT2;
dist(:,:,3) = im_pre_lef + im_nex_rig - BT2; 
dist(:,:,4) = im_nex_lef + im_pre_rig - BT2;
dist(:,:,5) = tmp1 + im_r_lef + im_pre_lef; 
dist(:,:,6) = tmp1 + im_r_rig + im_pre_rig;
dist(:,:,7) = tmp2 + im_r_lef + im_nex_lef; 
dist(:,:,8) = tmp2 + im_r_rig + im_nex_rig;
dist(:,:,1) = 1.5*dist(:,:,1); %% scale to the same level
dist(:,:,2) = 1.5*dist(:,:,2);
dist(:,:,3) = 1.5*dist(:,:,3);
dist(:,:,4) = 1.5*dist(:,:,4);
%find the signed distance with minimal absolute value
tmp = abs(dist); 
[~,ind] = min(tmp,[],3); 
%turn sub to index, but faster than sub2ind
dim1 = uint32(size(dist,1));
dim2 = uint32(size(dist,1)*size(dist,2));
index = row + dim1.*uint32(col-1) + dim2.*uint32(ind-1);
dm = single(step/3)*dist(index); 
%update current pixels
res(BT_r,BT_c) = BT + dm;

%% %%%%%%%%%%%%%%%%%%%% curvature energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%these energies are just for observation (so, formula can be changed by user)
function en = curv_TV(im)
[gx,gy]=mygrad(im); 
g = abs(gx) + abs(gy);
en = sum(g(:));
function en = curv_MC(im)
[gx,gy]=mygrad(im);
[gxx,gxy]=mygrad(gx);
[gyx,gyy]=mygrad(gy);
%standard scheme
num = (1+gy.^2).*gxx - gx.*gy.*(gxy+gyx)+ (1+gx.^2).*gyy;
den = (1+gx.^2+gy.^2);
den = sqrt(den).*den*2;
g = num./den;
en = sum(abs(g(:)));
function en = curv_GC(im)
[gx,gy]=mygrad(im);
[gxx,gxy]=mygrad(gx);
[gyx,gyy]=mygrad(gy);
%standard scheme
num = gxx.*gyy-gxy.*gyx;
den = 1+gx.^2+gy.^2;
den = den.*den;
g = num./den; 
en = sum(abs(g(:)));
function [gx, gy]=mygrad(im)
gx=[im(:,2)-im(:,1) (im(:,3:end)-im(:,1:end-2))./2 im(:,end)-im(:,end-1)];
gy=[im(2,:)-im(1,:) ; (im(3:end,:)-im(1:end-2,:))./2 ; im(end,:)-im(end-1,:)];