function [f] = get_figure(pixel,fig,image_size)
[x,y] = ind2sub(image_size,pixel);
a = repmat([x y], [size(fig,1) 1]) + fig;
ind = (a(:,1)>0 & a(:,1)<= image_size(1) & a(:,2)>0 & a(:,2)<= image_size(2));
a = a(ind,:);
f = sub2ind(image_size,a(:,1),a(:,2));
end