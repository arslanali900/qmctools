function result = hn(n, r)
  L = 3.0;
  rc = 0.5*L;
  omega = L*L*L;
  M = 5;
  delta = rc/(M-1);
  i = (n-1)/3;
  alpha = mod(n-1,3);
  ra = delta*(i-1);
  rb = delta*i;
  rc = delta*(i+1);
  S = [ 1.0   0.0   0.0  -10.0  15.0  -5.0; ...
        0.0   1.0   0.0   -6.0   8.0  -3.0; ...
        0.0   0.0   0.5   -1.5   1.5  -0.5 ];
  
  if ((r > rb) & (r <= rc)) 
    sum = 0;
    prod = 1;
    for j = [1:6] 
      sum = sum + S(alpha+1,j)*prod;
      prod = prod * ((rb-r)/delta);
    end;
    for j = [0:alpha-1]
      sum = sum * (-delta);
    end;
      result = sum;
  elseif ((r > rb) && (r <= rc))
    sum = 0;
    prod = 1.0;
    for j = [1:6]
      sum = sum + S(alpha+1,j)*prod;
      prod = prod * ((r-rb)/delta);
    end;
    for j = [0:alpha-1]
      sum = sum * delta;
    end;
    result = sum;
  else
    result = 0;
  end;
  
  
  
