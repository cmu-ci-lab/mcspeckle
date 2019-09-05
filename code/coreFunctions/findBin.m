function [bin_ind,bin_min,bin_max]=findBin(x,box_min,box_max,box_bin);


Nb=round((box_max-box_min)./box_bin);

ii=ceil((x-box_min)./box_bin);


bin_ind=ii(1);
for j=2:length(ii)
 bin_ind=bin_ind+(ii(j)-1)*prod(Nb(1:j-1));
end
if nargout>1
bin_min=box_min+box_bin.*(ii-1);
bin_max=bin_min+box_bin;


end