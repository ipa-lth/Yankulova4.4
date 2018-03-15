function plot_moving_poles(A, b, c, d, k_0, k_1, pmin)
  poles = [];
  for p = 1:-0.001:pmin
      sys_closed = ss(A-b*(k_0+(1.0/p)*k_1)', b, c, d);
      pol = pole(sys_closed);
      poles = [poles, pol];
  end
   
    # another approach to plot
    real_part = real(poles);
    imag_part = imag(poles);

    # Display a window with a plot of real, imag
    plot(real_part', imag_part', 'b-',
         real_part(:,1), imag_part(:,1), 'b*',
         real_part(:, size(real_part,2)), imag_part(:, size(imag_part,2)), 'rx');
return