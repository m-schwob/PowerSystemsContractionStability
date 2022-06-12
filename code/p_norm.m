function P_norm = p_norm (P,d,w,d_eq,w_eq)
    x = P*[(d-d_eq);(w-w_eq)];
    P_norm = norm(x);
end
