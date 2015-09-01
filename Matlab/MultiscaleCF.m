function r = MultiscaleCF(im, ScaleNum, FilterType, ItNum)
%%% im = imput image; ScaleNum = number of scale space; FilterTyep (1=MC, 2 = GC); 
%%% ItNum = iteration on each scale

pyramid = cell(ScaleNum,1);
pyramid{1} = double(im);
for i=2:ScaleNum
    [A,H,V,D]= lwt2(pyramid{i-1},'sym4');
    pyramid{i-1} = [A, H; V, D];
    pyramid{i} = A;
end

for i=ScaleNum:-1:2
    if (FilterType == 1) 
        A = MCF(pyramid{i}, ItNum);
        [m,n]=size(pyramid{i-1});
        H = pyramid{i-1}(1:m/2,n/2+1:end);
        V = pyramid{i-1}(m/2+1:end,1:n/2);
        D = pyramid{i-1}(m/2+1:end,n/2+1:end);
        pyramid{i-1} = ilwt2(A,H,V,D,'sym4');
    else
        A = GCF(pyramid{i}, ItNum);
         [m,n]=size(pyramid{i-1});
        H = pyramid{i-1}(1:m/2,n/2+1:end);
        V = pyramid{i-1}(m/2+1:end,1:n/2);
        D = pyramid{i-1}(m/2+1:end,n/2+1:end);
        pyramid{i-1} = ilwt2(A,H,V,D,'sym4');
    end
end

i = 1;
if (FilterType == 1) 
    pyramid{i} = MCF(pyramid{i}, ItNum);
else
    pyramid{i}  = GCF(pyramid{i}, ItNum);
end

r = pyramid{1};


%************* internal functions*****************%

function r = MCF(im, ItNum)
im = double(im); [m,n]=size(im);  r = im;
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
    r = MCUpdate(r,BC,BC_pre,BC_nex,BC_lef,BC_rig,BC_lu,BC_ld,BC_ru,BC_rd);
    r = MCUpdate(r,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd);
    r = MCUpdate(r,WC,WC_pre,WC_nex,WC_lef,WC_rig,WC_lu,WC_ld,WC_ru,WC_rd);
    r = MCUpdate(r,WT,WT_pre,WT_nex,WT_lef,WT_rig,WT_lu,WT_ld,WT_ru,WT_rd);
end
function res = MCUpdate(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd)
res = im; BT8 = 8*im(BT); 
dist = zeros(size(BT_pre,1),4);
tmp1 = 2.5*(im(BT_pre) + im(BT_nex)) - BT8;
tmp2 = 2.5*(im(BT_lef) + im(BT_rig)) - BT8;

dist(:,1) = tmp1  + 5*im(BT_rig) - im(BT_ru) - im(BT_rd);
dist(:,2) = tmp1  + 5*im(BT_lef) - im(BT_lu) - im(BT_ld);
dist(:,3) = tmp2  + 5*im(BT_pre) - im(BT_lu) - im(BT_ru);
dist(:,4) = tmp2  + 5*im(BT_nex) - im(BT_ld) - im(BT_rd);

dist(:,1:4) = dist(:,1:4)/8; %% minimal projection
dist= dist'; tmp = abs(dist); [v,ind] = min(tmp);
tmp = sub2ind(size(dist),ind',(1:size(dist,2))');
res(BT) = res(BT) + dist(tmp);

function r = GCF(im, ItNum)
im = double(im); [m,n]=size(im); r = im;
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
    r = GCUpdate(r,BC,BC_pre,BC_nex,BC_lef,BC_rig,BC_lu,BC_ld,BC_ru,BC_rd);
    r = GCUpdate(r,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd);
    r = GCUpdate(r,WC,WC_pre,WC_nex,WC_lef,WC_rig,WC_lu,WC_ld,WC_ru,WC_rd);
    r = GCUpdate(r,WT,WT_pre,WT_nex,WT_lef,WT_rig,WT_lu,WT_ld,WT_ru,WT_rd);
end
function res = GCUpdate(im,BT,BT_pre,BT_nex,BT_lef,BT_rig,BT_lu,BT_ld,BT_ru,BT_rd)
res = im; BT2 = 2*im(BT); BT3 = 3*im(BT);
dist = zeros(size(BT_pre,1),8);
dist(:,1) = im(BT_pre) + im(BT_nex) - BT2; dist(:,2) = im(BT_lef) + im(BT_rig) - BT2;
dist(:,3) = im(BT_lu) + im(BT_rd) - BT2; dist(:,4) = im(BT_ld) + im(BT_ru) - BT2;

dist(:,5) = im(BT_pre) + im(BT_lef) + im(BT_lu) - BT3; dist(:,6) = im(BT_pre) + im(BT_rig) + im(BT_ru) - BT3;
dist(:,7) = im(BT_nex) + im(BT_lef) + im(BT_ld) - BT3; dist(:,8) = im(BT_nex) + im(BT_rig) + im(BT_rd) - BT3;

dist(:,1:4) = dist(:,1:4)/2; dist(:,5:8) = dist(:,5:8)/3; %% minimal projection
dist= dist'; tmp = abs(dist); [v,ind] = min(tmp);
tmp = sub2ind(size(dist),ind',(1:size(dist,2))');
res(BT) = res(BT) + dist(tmp);
