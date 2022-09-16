% gives the matrix measure of any point of mesh grid (3d)
function MU = grid_matmis3(J_in, delta_2, omega_1, omega_2,L)
len_x = length(omega_1);
len_y = length(omega_2);
len_z = length(delta_2);
MU = zeros (len_x,len_y,len_z);
    for w1 = 1:len_x
        for w2 =1:len_y
            for d2 = 1:len_z
                MU(w1,w2,d2) = matmis(J_in(delta_2(w1,w2,d2), omega_1(w1,w2,d2), omega_2(w1,w2,d2)), L);
            end
        end
    end
end
