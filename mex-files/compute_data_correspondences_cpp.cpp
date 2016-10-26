//mex mex-files\compute_correspondences_cpp.cpp -IE:\OneDrive\EPFL\Code\external\eigen_dir
#define AU_MEX_UNCHECKED
#include "au_mex.h"
#include <Eigen/Dense>

Eigen::Vector2f projection_segment(Eigen::Vector2f p, Eigen::Vector2f c1, Eigen::Vector2f c2) {
    Eigen::Vector2f s;
    
    Eigen::Vector2f u = c2 - c1;
    Eigen::Vector2f v = p - c1;
    
    float alpha = u.dot(v) / u.dot(u);
    
    if (alpha <= 0)
        s = c1;
    
    if (alpha >= 1)
        s = c2;
    
    if (alpha > 0 && alpha < 1)
        s = c1 + alpha * u;
    return s;
}

void compute_correspondences(mwSize num_data, mwSize num_blocks,
        const mlx_array<mlx_double> & C, const mlx_array<mlx_double> & B, const mlx_array<mlx_double> & D,
        mlx_array<mlx_double> & M, mlx_array<mlx_double> & I) {
    
    for(int i = 0; i < num_data; i++) {
        Eigen::Vector2f p = Eigen::Vector2f(D(i, 0), D(i, 1));
        float min_distance = RAND_MAX;
        
        for(int j = 0; j < num_blocks; j++) {
            
            mwIndex index1 = (mwIndex) B(j, 0) - 1;
            mwIndex index2 = (mwIndex) B(j, 1) - 1;
            Eigen::Vector2f c1 = Eigen::Vector2f(C(index1, 0), C(index1, 1));
            Eigen::Vector2f c2 = Eigen::Vector2f(C(index2, 0), C(index2, 1));
            
            Eigen::Vector2f q = projection_segment(p, c1, c2);
            
            float distance = (p - q).norm();
            if (distance < min_distance) {
                min_distance = distance;
                M(i, 0) = q[0]; M(i, 1) = q[1];
                I(i, 0) = j + 1;
            }
        }
    }
}

void mlx_function(mlx_inputs& in, mlx_outputs& out) {

    mlx_array<mlx_double> sizes_array(mlx_size{1, 3}, in[0]);
    mwSize num_centers = (mwSize) sizes_array[0];
    mwSize num_blocks = (mwSize) sizes_array[1];
    mwSize num_data = (mwSize) sizes_array[2];
    
    mlx_array<mlx_double> C(mlx_size{num_centers, 2}, in[1]);
    mlx_array<mlx_double> B(mlx_size{num_blocks, 2}, in[2]);
    mlx_array<mlx_double> D(mlx_size{num_data, 2}, in[3]);
    
    mlx_make_array<double> M(mlx_size{num_data, 2});
    mlx_make_array<double> I(mlx_size{num_data, 1});
    
    compute_correspondences(num_data, num_blocks, C, B, D, M, I);
    
    out[0] = I;
    out[1] = M;
}

