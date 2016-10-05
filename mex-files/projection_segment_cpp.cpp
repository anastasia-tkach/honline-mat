//mex mex-files\projection_segment_cpp.cpp -IE:\OneDrive\EPFL\Code\external\eigen_dir
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

void mlx_function(mlx_inputs& in, mlx_outputs& out)
{
    mlx_array<mlx_double> p_array(mlx_size{2,1}, in[0]);
    mlx_array<mlx_double> c1_array(mlx_size{2,1}, in[1]);
    mlx_array<mlx_double> c2_array(mlx_size{2,1}, in[2]);
    mlx_make_array<double> s_array(p_array.size);
    
    // Parse inputs to eigen
    Eigen::Vector2f p = Eigen::Vector2f(p_array[0], p_array[1]);
    Eigen::Vector2f c1 = Eigen::Vector2f(c1_array[0], c1_array[1]);
    Eigen::Vector2f c2 = Eigen::Vector2f(c2_array[0], c2_array[1]);
    
    // Do the math
    Eigen::Vector2f s = projection_segment(p, c1, c2);
       
    // Parse outputs from eigen
    for(mwSize i = 0; i < s_array.numel(); ++i) 
        s_array[i] = s(i);    
    
    out[0] = s_array;
}

