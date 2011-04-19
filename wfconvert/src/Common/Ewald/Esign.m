function result = Esign (sign, delta, omega, j, k, n)
  if (n == 0)
    result = -sign * 4*pi*i/(omega*k*k)*(exp(sign*i*k*delta)-1.0)* ...
             exp(i*k*j*delta);
  else
    result = -i/k*(4*pi/(omega*k)*sign^(n+1)*exp(i*k*delta*(j+sign)) ...
                   - n/delta * Esign(sign, delta, omega, j, k, n- ...
                                     1));
  end;
  
    
