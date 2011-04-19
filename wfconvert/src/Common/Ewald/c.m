function result = c(n, k)
  result = zeros(size(k));
  for j = [1:length(k)]
    result(j) = cone(n,k(j));
  end;
  
