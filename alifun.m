function [ma mi] = alifun(A)
  aw = A(:);
  ma = aw(1);
  mi = aw(1);
  for  i1 = 2:numel(aw)
     if ma < aw(i1), ma = aw(i1);
     elseif mi > aw(i1), mi = aw(i1);  
     end
  end
end