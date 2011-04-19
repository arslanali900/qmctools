%bands = load ('TestBands.dat');
bands = load ('TestLDABands32k_kcut10.dat');
bands = bands(:,1:9);
bands2 = load ('TestBands2.dat');
t = [0:0.1:6.9];
t2 = [0:0.1:6.9];

a = 0.5*8.011*[-1 1 1; 1 -1 1; 1 1 -1];
b = 2*pi/8.011*[0 1 1; 1 0 1; 1 1 0];


bands0  = bands(1,4);
bands02  = bands2(1,4);

G1 = ones(70,1)*[0      2*pi/8     0      ];
G2 = ones(70,1)*[0      0          2*pi/8 ];
G3 = ones(70,1)*[2*pi/8 0          0      ];
G4 = ones(70,1)*[2*pi/8 2*pi/8     0      ];
G5 = ones(70,1)*[0      2*pi/8     2*pi/8 ];
G5 = ones(70,1)*[0      2*pi/8     2*pi/8 ];
%G1 = ones(70,1)*[2*pi/8.011 0 0];

bands(:,4:9) = 2*(bands(:,4:9)-bands0);
bands2(:,4:9) = 2*(bands2(:,4:9)-bands02);
fp1 = sum(bands(:,1:3).^2,2);
fp2 = sum((bands(:,1:3)-G1).^2,2);
fp3 = sum((bands(:,1:3)-G2).^2,2);
fp4 = sum((bands(:,1:3)-G3).^2,2);
fp5 = sum((bands(:,1:3)-G4).^2,2);
fp6 = sum((bands(:,1:3)-G5).^2,2);


hold off;

abands = load ('ABINITBands.dat');
tb = [0:1:49]/49*7;
abands0 = abands(1,4);
abands(:,4:11) = 2*(abands(:,4:11)-abands0);


handles = plot (t, bands(:,4), 'b', tb, abands(:,4), 'r', t, fp1, 'k', ...% t2, bands2(:,4), 'g--',...
      t, bands(:,5), 'b', tb, abands(:,5), 'r', t, fp2, 'k', ...% t2, bands2(:,5), 'g--',...
      t, bands(:,6), 'b', tb, abands(:,6), 'r', t, fp3, 'k', ...% t2, bands2(:,6), 'g--',...
      t, bands(:,7), 'b', tb, abands(:,7), 'r', t, fp4, 'k', ...% t2, bands2(:,7), 'g--',...
      t, bands(:,8), 'b', tb, abands(:,8), 'r', t, fp5, 'k', ...% t2, bands2(:,8), 'g--',...
      t, bands(:,9), 'b', tb, abands(:,9), 'r', t, fp6, 'k');
set (handles, 'LineWidth', 2);
axis([0 7 0 0.55]);

ca = get (gcf, 'CurrentAxes');
set(ca, 'FontSize', 16');
xlabel ('k Points');
ylabel ('Energy (Ry)');


set(ca, 'XTick', [0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0]);
set(ca, 'XTickLabel', ['G';'H';'N';'G';'P';'N';'P';'H']);
legend ('PH bands', 'ABINIT bands', 'Free particle bands');
%title ('BCC Sodium Band Comparison');
