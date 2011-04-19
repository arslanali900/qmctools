function result = cone(n, k)
  L = 8.0;
  rc = 0.5*L;
  omega = L*L*L;
  M = 5;
  delta = rc/(M-1);
  i = floor((n-1)/3);
  alpha = mod(n-1,3);
  S = [ 1.0,   0.0,   0.0,  -10.0,  15.0,  -6.0; ...
        0.0,   1.0,   0.0,   -6.0,   8.0,  -3.0; ...
        0.0,   0.0,   0.5,   -1.5,   1.5,  -0.5 ];
  sum = 0.0;
  for n=[0:5]
    sum = sum + S(alpha+1, n+1)*...
          (Dsign(1,  delta, omega, i, k, n) +...
           (-1)^(alpha+n)*Dsign(-1, delta, omega, i, k, n));
  end;
  result = delta^alpha * sum;
  
