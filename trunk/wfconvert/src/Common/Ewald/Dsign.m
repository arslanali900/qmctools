function result = Dsign (sign, delta, omega, j, k, n)
  result = delta * imag(Esign(sign, delta, omega, j, k, n+1) + ...
                        i*Esign(sign, delta, omega, j, k, n));
  
