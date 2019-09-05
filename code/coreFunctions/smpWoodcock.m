function d=smpWoodcock(x,box_min,box_max,box_bin,sigt,v)
      sigt_max=max(sigt(:));
      
      d=0;
      while 1
          d=d-log(1-rand)/sigt_max;
          nx=x+d*v;
          if(max(nx>box_max)||max(nx<box_min))
              break
          end
          [bin_ind]=findBin(nx,box_min,box_max,box_bin);
          if rand<=(sigt(bin_ind)/sigt_max)
              break
          end
      end
end