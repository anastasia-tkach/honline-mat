//mex jacobian_shape_pose_cpp.cpp -IE:\OneDrive\EPFL\Code\external\eigen_dir
#define AU_MEX_UNCHECKED
#include "au_mex.h"
#include <Eigen/Dense>


void jacobian_pose(const mlx_array<mlx_double> & beta, const mlx_array<mlx_double> & theta,
        mwSize num_points, mwSize num_segments, mwSize num_joints, mwSize max_kinematic_chain, const mlx_array<mlx_double> & DataPoints,
        const mlx_array<mlx_double> & ModelPoints, const mlx_array<mlx_double> & segment_indices, const mlx_array<mlx_double> & SegmentsKinematicChain,
        const mlx_array<mlx_double> & SegmentsGlobal, const mlx_array<mlx_double> & JointsSegmentId, const mlx_array<mlx_double> & JointsAxis,
        mlx_array<mlx_double> & F, mlx_array<mlx_double> & J) {
    
    for (int k = 0; k < num_points; k++) {
        //mexPrintf("\nk = %d\n", k);
        Eigen::Vector3d d = Eigen::Vector3d(DataPoints(k, 0), DataPoints(k, 1), 0);
        Eigen::Vector3d m = Eigen::Vector3d(ModelPoints(k, 0), ModelPoints(k, 1), 0);
        Eigen::Vector3d n = (d - m) / (d - m).norm();
        
        Eigen::Matrix<double, 3, Eigen::Dynamic> j = Eigen::Matrix<double, 3, Eigen::Dynamic> ::Zero(3, num_segments - 1 + num_joints - 1);
        
        for (int l = 0; l < max_kinematic_chain; l++) {            
            if (SegmentsKinematicChain(segment_indices(k, 0) - 1, l) == -1) break;
            
            //mexPrintf("l = %d\n", l);
            
            mwIndex joint_id = (mwIndex) SegmentsKinematicChain(segment_indices(k, 0) - 1, l) - 1;
            mwIndex segment_id = (mwIndex) JointsSegmentId(joint_id, 0) - 1;
            Eigen::Vector4d u = Eigen::Vector4d(JointsAxis(joint_id, 0), JointsAxis(joint_id, 1), JointsAxis(joint_id, 2), 1);
            
            Eigen::Vector3d p = Eigen::Vector3d(SegmentsGlobal(segment_id, 12), SegmentsGlobal(segment_id, 13), SegmentsGlobal(segment_id, 14));
            
            Eigen::Matrix<double, 4, 4> T = Eigen::Matrix<double, 4, 4>();
            for (int x = 0; x < 4; x++) {
                for (int y = 0; y < 4; y++) {
                    T(x, y) = SegmentsGlobal(segment_id, 4 * y + x);
                    //mexPrintf("%f  ", T(i, j));
                }
                //mexPrintf("\n");
            }
            //mexPrintf("\n");
            
            // shape
            Eigen::Vector4d w4  = T * Eigen::Vector4d(0, 1, 0, 1);
            Eigen::Vector3d w = Eigen::Vector3d(w4(0) / w4(3), w4(1) / w4(3), w4(2) / w4(3));
            w = w - p;
            
            double c = 1;
            bool last_in_kinematic_chain = SegmentsKinematicChain(segment_indices(k, 0) - 1, l + 1) == -1;
            if (l == 2 || last_in_kinematic_chain)
                c = (m - p).norm() / beta(l, 0);
            j.col(segment_id) = c * w;
            
            // pose
            Eigen::Vector4d v4 = T * u;
            Eigen::Vector3d v = Eigen::Vector3d(v4(0) / v4(3), v4(1) / v4(3), v4(2) / v4(3));
            v = v - p;
            j.col(num_segments - 1 + joint_id) = v.cross(m - p);
            
            
            /*mexPrintf("w = %f %f %f\n", w(0), w(1), w(2));
            mexPrintf("c = %f\n", c);
            mexPrintf("v = %f %f %f\n", v(0), v(1), v(2));
            mexPrintf("m = %f %f %f\n", m(0), m(1), m(2));
            mexPrintf("p = %f %f %f\n", p(0), p(1), p(2));
            mexPrintf("n = %f %f %f\n", n(0), n(1), n(2));
            mexPrintf("d = %f %f %f\n", d(0), d(1), d(2));*/
            
        }
        /*mexPrintf("j = \n");
        for (int y = 0; y < 3; y++) {
            for (int x = 0; x < num_segments - 1 + num_joints - 1; x++) {             
                mexPrintf("%f\t", j(y, x));
            }
            mexPrintf("\n ");
        }*/        
        
        for (int var = 0; var < num_segments - 1 + num_joints - 1; var++) {
            J(k, var) = - n.dot(j.col(var));
            //mexPrintf("%f\t", J(k, var));
        }
        //mexPrintf("\n");
        F(k, 0) = n.dot(d - m);
    }
    
}

void mlx_function(mlx_inputs& in, mlx_outputs& out) {
    //(sizes, DataPoints, ModelPoints, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis)
    
    mlx_array<mlx_double> beta(mlx_size{3, 1}, in[0]);
    mlx_array<mlx_double> theta(mlx_size{3, 1}, in[1]);
    
    mlx_array<mlx_double> sizes_array(mlx_size{1, 4}, in[2]);
    mwSize num_points = (mwSize) sizes_array[0];
    mwSize num_joints = (mwSize) sizes_array[1];
    mwSize num_segments = (mwSize) sizes_array[2];
    mwSize max_kinematic_chain = (mwSize) sizes_array[3];
    
    mlx_array<mlx_double> DataPoints(mlx_size{num_points, 2}, in[3]);
    mlx_array<mlx_double> ModelPoints(mlx_size{num_points, 2}, in[4]);
    mlx_array<mlx_double> segment_indices(mlx_size{num_points, 1}, in[5]);
    mlx_array<mlx_double> SegmentsKinematicChain(mlx_size{num_segments, max_kinematic_chain}, in[6]);
    mlx_array<mlx_double> SegmentsGlobal(mlx_size{num_segments, 16}, in[7]);
    mlx_array<mlx_double> JointsSegmentId(mlx_size{num_joints, 1}, in[8]);
    mlx_array<mlx_double> JointsAxis(mlx_size{num_joints, 3}, in[9]);
    
    mlx_make_array<double> F(mlx_size{num_points, 1});
    mlx_make_array<double> J(mlx_size{num_points, num_segments - 1 + num_joints - 1});
    
    /*mlx_make_array<double> DataPoints_(mlx_size{num_points, 2});
     * mlx_make_array<double> ModelPoints_(mlx_size{num_points, 2});
     * mlx_make_array<double> segment_indices_(mlx_size{num_points, 1});
     * mlx_make_array<double> SegmentsKinematicChain_(mlx_size{num_segments, max_kinematic_chain});
     * mlx_make_array<double> SegmentsGlobal_(mlx_size{num_segments, 16});
     * mlx_make_array<double> JointsSegmentId_(mlx_size{num_joints, 1});
     * mlx_make_array<double> JointsAxis_(mlx_size{num_joints, 3});
     * out[0] = DataPoints_;
     * out[1] = ModelPoints_;
     * out[2] = segment_indices_;
     * out[3] = SegmentsKinematicChain_;
     * out[4] = SegmentsGlobal_;
     * out[5] = JointsSegmentId_;
     * out[6] = JointsAxis_;*/
    
    jacobian_pose(beta, theta, num_points, num_segments, num_joints, max_kinematic_chain, DataPoints, ModelPoints, segment_indices,
            SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis, F, J);
    
    out[0] = F;
    out[1] = J;
}

