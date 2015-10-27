function [result, Energy] = TVFilter(im, ItNum)
im = single(im); [m,n]=size(im); Energy = zeros(ItNum,1); result = im;
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
%% dual Mesh optimization
for i = 1:ItNum
    Energy(i) = TVEnergy(result);
    result = SimpleUpdate(result,BC,BC_pre,BC_nex,BC_lef,BC_rig,BC_lu,BC_ld,BC_ru,BC_rd);
    result = SimpleUpdate(result,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd);
    result = SimpleUpdate(result,WC,WC_pre,WC_nex,WC_lef,WC_rig,WC_lu,WC_ld,WC_ru,WC_rd);
    result = SimpleUpdate(result,WT,WT_pre,WT_nex,WT_lef,WT_rig,WT_lu,WT_ld,WT_ru,WT_rd);
end
function res = SimpleUpdate(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd)
res = im; BT5 = 5*im(BT); dist = zeros(size(BT_pre,1),8,'single');
dist(:,1) = im(BT_pre) + im(BT_nex) + im(BT_lu) + im(BT_lef) + im(BT_ld) - BT5; 
dist(:,2) = im(BT_pre) + im(BT_nex) + im(BT_ru) + im(BT_rig) + im(BT_rd) - BT5;
dist(:,3) = im(BT_lef) + im(BT_rig) + im(BT_lu) + im(BT_pre) + im(BT_ru) - BT5; 
dist(:,4) = im(BT_lef) + im(BT_rig) + im(BT_ld) + im(BT_nex) + im(BT_rd) - BT5; 
dist(:,5) = im(BT_pre) + im(BT_lef) + im(BT_lu) + im(BT_ru) + im(BT_ld) - BT5; 
dist(:,6) = im(BT_pre) + im(BT_rig) + im(BT_ru) + im(BT_lu) + im(BT_rd) - BT5;
dist(:,7) = im(BT_nex) + im(BT_lef) + im(BT_ld) + im(BT_lu) + im(BT_rd) - BT5; 
dist(:,8) = im(BT_nex) + im(BT_rig) + im(BT_rd) + im(BT_ld) + im(BT_ru) - BT5;

dist = dist/5; %% minimal projection
dist= dist'; tmp = abs(dist); [v,ind] = min(tmp);
tmp = sub2ind(size(dist),ind',(1:size(dist,2))');
res(BT) = res(BT) + dist(tmp);
function en = TVEnergy(im) %compute TV energy
t = double(im);[gx,gy]=gradient(t); g = abs(gx) + abs(gy);
en = sum(g(:));