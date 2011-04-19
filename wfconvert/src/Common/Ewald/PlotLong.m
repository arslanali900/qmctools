V = load ('Vlong.dat');
r = V(:,1);



alpha = 7.0/10;

Vewald = -r.^-1.*erf(alpha*r);
Vcoul = -r.^-1;
plot (V(:,1), V(:,2), r, Vewald); a = axis;
plot (r(1:5000), V(1:5000,2), r, Vewald, r, Vcoul);
axis(a);

xlabel ('r');
ylabel ('V(r)');
legend ('V_{optimized}(r)', 'V_{ewald}(r)', 'V_{bare}(r)',2);
title ('Optimized long range potential');
text (8, -0.75, date);
