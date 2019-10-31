function mask = sampling_mask(acceleration,resx,resy,pa,pb)

%place=2;%1 Mac desktop omega,2 Mac stick 3 PC work 4 PC home
iterations=10;

%acceleration0 =8
%acceleration =acceleration0*resxy/256

dims_xy=[resx resy];

paf=pa*10;
% pbtest=resxy*0.05;
%pb0= 5.6
%pb=pb0*resxy/256
%pb=~fraction * 0.05
pbf=pb*10;
if acceleration >= 1
sampling_fraction = 1/acceleration;
else
sampling_fraction = acceleration;
acceleration = 1/sampling_fraction;
end
if numel(dims_xy) == 1
dims_xy = [dims_xy dims_xy];
end

[mypdf,val] = genPDF_wn_v2(dims_xy,pa,sampling_fraction ,pb,1);

[mask, stat, actpctg]= genSampling(mypdf,iterations,1);
% mask=reshape(mask,resxy*resxy,1);
% mask=transpose(mask);
if ~mask(end)
    mask(find(mask,1,'last'))=0;
    mask(end)=1; 
end
end