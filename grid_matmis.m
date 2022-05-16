% gives the matrix measure of any point of mesh greid
function MU = grid_matmis(func_A, delta_2, omega_2, L)
    [x_size,y_size] = size(delta_2); 
    MU = zeros(x_size, y_size);
    for x = 1:x_size
        for y =1:y_size
            MU(x,y) = matmis(func_A(delta_2(x,y), omega_2(x,y)), L);
        end
    end
end
