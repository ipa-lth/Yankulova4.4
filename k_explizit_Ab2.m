function [k0, k1] = k_explizit_Ab2(roots_p1, roots_pmin, pmin, A, b)

    r0 = place(A, b, roots_p1); #k(p=1)
    r1 = place(A, b, roots_pmin); #k(p=pmin)

    # This seems to work as expected
    if(pmin < 1.0)
      k1 = 1.0/(1.0-1.0/pmin) * (r0 - r1);
    else
      k1 = r0 - r1;
    end
    k0 = r0 - k1;

    k0 = k0';
    k1 = k1';
    
return