function [Hind, E, I] = analyze_fft(Hz_ind0, Ifft0, Efft0)

Hz_bound = [1 100];

sample_ind = Hz_bound(1)<=Hz_ind0 & Hz_ind0<=Hz_bound(2);

Hz_ind = Hz_ind0(sample_ind);
Ifft1 = Ifft0(sample_ind);
Efft1 = Efft0(sample_ind);

% filter
gsig = 10;       % [Hz]
gfil = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
gfil = gfil/sum(gfil(:));

Efft2 = conv(Efft1, gfil, 'same');
Ifft2 = conv(Ifft1, gfil, 'same');

% figure
% hold on
% plot(Hz_ind, Ifft1, 'b-')
% plot(Hz_ind, Efft1, 'r-')
% plot(Hz_ind, Ifft2, 'b-', 'linewidth', 4)
% plot(Hz_ind, Efft2, 'r-', 'linewidth', 4)

Hind = Hz_ind(:);
E.FF = Efft2(:);
I.FF = Ifft2(:);

E.FFmaxHz = Hind(Efft2(:) == max(Efft2(:)));
I.FFmaxHz = Hind(Ifft2(:) == max(Ifft2(:)));

end
