function [pdf,val] = genPDF_wn_v2(imSize,pa,pctg,pb,disp)

%[pdf,val] = genPDF(imSize,p,pctg [,distType,radius,disp])
%
%	generates a pdf for a 1d or 2d random sampling pattern
%	with polynomial variable density sampling
%
%	Input:
%		imSize - size of matrix or vector
%		p - power of polynomial
%		pctg - partial sampling factor e.g. 0.5 for half
%		distType - 1 or 2 for L1 or L2 distance measure
%		radius - radius of fully sampled center
%		disp - display output
%
%	Output:
%		pdf - the pdf
%		val - min sampling density
%
% 
%	Example:
%	[pdf,val] = genPDF([128,128],2,0.5,2,0,1);
%
%	(c) Michael Lustig 2007

% imSize=[256 256]
% p=14
% pctg=0.125
% distType=2
% radius=0
% disp=1
% 
% 
minval=0;
maxval=1;
val = 0.5;
% 
% if length(imSize)==1
% 	imSize = [imSize,1];
% end
% 
% sx = imSize(1);
% sy = imSize(2);
% PCTG = floor(pctg*sx*sy);
% 
% % 
% % a=3;
% % b=2.5;
% 
% 
% if sum(imSize==1)==0  % 2D
% 	[x,y] = meshgrid(linspace(-1,1,sy),linspace(-1,1,sx));
% 	switch distType
% 		case 1
% 			r = max(abs(x),abs(y));
% 		otherwise
% 			r = ((sqrt(x.^2+y.^2)).^2.1);
% 			r = r/max(abs(r(:)));			
% 	end
% 
% else %1d
% 	r = abs(linspace(-1,1,max(sx,sy)));
% end
% figure;imshow(r,[])
% 
% 
% 
% idx = find(r<radius);
% 
% pdf = (1-r).^p;
% pdf(idx) = 1;
sx = imSize(1);
sy = imSize(2);
PCTG = floor(pctg*sx*sy);
iraw=zeros(imSize);
xres=size(iraw,1); yres=size(iraw,2); zres=size(iraw,3);
[y,x,z] = meshgrid(-yres/2:yres/2-1,-xres/2:xres/2-1,-zres/2:zres/2-1);

x2=xres/2; y2=yres/2; z2=zres/2;

% f1=exp(-((pb*sqrt(x.^2)/sx).^pa));
% f2=exp(-((pb*sqrt(y.^2)/sy).^pa));

% f1=sqrt(exp(-((pb*sqrt(x.^2)/sx).^pa)));
% f2=sqrt(exp(-((pb*sqrt(y.^2)/sy).^pa)));

f1=sqrt(exp(-((pb*sqrt(x.^2+y.^2)/sx).^pa)));
f2=sqrt(exp(-((pb*sqrt(y.^2+x.^2)/sy).^pa)));


f=f1.*f2; 
f=f/max(f(:));
% figure;imshow(f,[])
% figure;plot(1:sx,f1)
pdf=f;
if floor(sum(pdf(:))) > PCTG
	error('infeasible without undersampling dc, change pa or pb');
end
% 
% % begin bisection
while(1)
	val = minval/2 + maxval/2;
	pdf = f + val; pdf(find(pdf>1)) = 1;
	N = floor(sum(pdf(:)));
	if N > PCTG % infeasible
		maxval=val;
	end
	if N < PCTG % feasible, but not optimal
		minval=val;
	end
	if N==PCTG % optimal
		break;
	end
end





if disp
	figure, 
	subplot(211), imshow(pdf)
	if sum(imSize==1)==0
		subplot(212), plot(pdf(end/2+1,:));
	else
		subplot(212), plot(pdf);
	end
end


% [mask,stat,actpctg] = genSampling(pdf,10,2);
%  size(find(mask==1))




