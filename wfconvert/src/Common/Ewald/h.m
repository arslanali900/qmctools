function result = h(n, r)
  result = zeros(size(r));
  for i=[1:length(r)]
    result(i) = hone(n, r(i));
  end;
  
  
  
