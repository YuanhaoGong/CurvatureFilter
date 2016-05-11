function res = HalfWindow(input, r, sigma, ItNum, stepsize)
%%% half window Gaussian filter
switch nargin
    case 1
        r = 2; sigma = 2; ItNum = 1; stepsize = 1;
    case 2
        sigma = 3; ItNum = 1; stepsize = 1;
    case 3
        ItNum = 1; stepsize = 1;
    case 4
        stepsize = 1;
    case 5
   otherwise
      disp('input is not correct.'), return;
end

%four half kernels by three separable kernel
h = fspecial('gaussian', 2*r+1, sigma); h = sum(h); 
h_v=h;h_v(r+2:end)=0; h_v = h_v/sum(h_v);
h_vf=fliplr(h_v); 
  
input = single(input); res = input;
for i=1:size(input,3)
    im = input(:,:,i);
    %pad the image
    im = padarray(im,[r,r],'replicate'); im = single(im); [m,n]=size(im);
    [col,row]=meshgrid(1:n,1:m); d = zeros(m,n,4,'single'); 
    ch = im;
    for it = 1:ItNum
        d(:,:,1) = conv2(h, h_v, ch,'same') - ch; 
        d(:,:,2) = conv2(h, h_vf, ch,'same') - ch;
        d(:,:,3) = conv2(h_v, h, ch,'same') - ch; 
        d(:,:,4) = conv2(h_vf, h, ch,'same') - ch;
        %find the minimal signed distance
        tmp = abs(d); 
        [v,ind] = min(tmp,[],3); 
        index = row + m*(col-1)+m*n*(ind-1);
        dm = d(index);
        %update
        ch = ch + stepsize*dm; 
    end
    %unpad
    res(:,:,i) = ch(r+1:end-r,r+1:end-r);
end
