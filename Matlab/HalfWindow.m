function res = HalfWindow(input, ItNum, step, r, sigma)
%perform half window Gaussian filter
input = single(input); res = input;
for i=1:size(input,3)
    im = input(:,:,i);
    im = padarray(im,[r,r],'replicate'); im = single(im); [m,n]=size(im);  ch = im;
    [col,row]=meshgrid(1:n,1:m); row = reshape(row,m*n,1); col = reshape(col,m*n,1); d = zeros(m,n,4,'single'); 
    h = fspecial('gaussian', 2*r+1, sigma); h = sum(h); h_v=h;h_v(r+2:end)=0; h_v = h_v/sum(h_v);h_vf=fliplr(h_v); 
    for it = 1:ItNum
        d(:,:,1) = conv2(h, h_v, ch,'same') - ch; d(:,:,2) = conv2(h, h_vf, ch,'same') - ch;
        d(:,:,3) = conv2(h_v, h, ch,'same') - ch; d(:,:,4) = conv2(h_vf, h, ch,'same') - ch;
        tmp = abs(d); [v,ind] = min(tmp,[],3); ind2 = reshape(ind,m*n,1);ind3 = sub2ind([m,n,4],row,col,ind2); 
        dm = d(ind3); dm = reshape(dm,m,n); ch = ch + step*dm; 
    end
    res(:,:,i) = ch(r+1:end-r,r+1:end-r);
end
