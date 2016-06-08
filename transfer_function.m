function output = transfer_function(output,eq_it)

% extract the base and surface motion

bs = output.defualtnodes(end); % base node

F1 = sprintf('%s%s%s','acc_sf = output.iteration.it_',num2str(eq_it),'.node_1.relative.TDVA(:,[1 4]);');
F2 = sprintf('%s%s%s%s%s','acc_bs = output.iteration.it_', num2str(eq_it),'.node_',num2str(bs),'.relative.TDVA(:,[1 4]);');

eval(F1);
eval(F2);


% taking the fft

dt = acc_sf(2,1)-acc_sf(1,1);
FS = 1/dt;
N  = size(acc_sf,1);
NP = 2^(nextpow2(N));
f  = (0:1/(NP-1):1)*FS;

fft_sf = fft(acc_sf(:,2),NP)/(N);
fft_bs = fft(acc_bs(:,2),NP)/(N);


transfer_function = abs(fft_sf)./abs(fft_bs);

TF = [f' transfer_function];


F3 = sprintf('%s%s%s','output.iteration.it_', num2str(eq_it), '.TF = TF;');
eval(F3)


end