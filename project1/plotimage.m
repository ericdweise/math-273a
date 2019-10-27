function [A] = plotimage(A)

   image(A'*255+1)
   Map = [0:1/255:1]'*ones(1,3);
   colormap(Map)
   axis('xy')
   axis('image')

end
