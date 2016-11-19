//mex jacobian_cpp.cpp -IE:\OneDrive\EPFL\Code\external\eigen_dir -IC:\Developer\include\glm
#define AU_MEX_UNCHECKED
#include "au_mex.h"
#include <Eigen/Dense>
#include "glm.hpp"


void jacobian_pose(const mlx_array<mlx_double> & beta, const mlx_array<mlx_double> & theta,
        mwSize num_points, mwSize num_segments, mwSize num_joints, mwSize max_kinematic_chain, const mlx_array<mlx_double> & DataPoints,
        const mlx_array<mlx_double> & ModelPoints, const mlx_array<mlx_double> & segment_indices, const mlx_array<mlx_double> & SegmentsKinematicChain,
        const mlx_array<mlx_double> & SegmentsGlobal, const mlx_array<mlx_double> & JointsSegmentId, const mlx_array<mlx_double> & JointsAxis,
        mlx_array<mlx_double> & F, mlx_array<mlx_double> & J) {
    
    for (int k = 0; k < num_points; k++) {
    //for (int k = 0; k < 1; k++) {
        //mexPrintf("\nk = %d\n", k);
        glm::dvec3 d = glm::dvec3(DataPoints(k, 0), DataPoints(k, 1), DataPoints(k, 2));
        glm::dvec3 m = glm::dvec3(ModelPoints(k, 0), ModelPoints(k, 1), ModelPoints(k, 2));
        glm::dvec3 n = (d - m) / glm::length(d - m);
        
        //mexPrintf("d = %f %f %f\n", d(0), d(1), d(2));
        //mexPrintf("m = %f %f %f\n", m(0), m(1), m(2));
        
        Eigen::Matrix<double, 3, Eigen::Dynamic> j = Eigen::Matrix<double, 3, Eigen::Dynamic> ::Zero(3, num_segments - 1 + num_joints - 1);
        
        for (int l = 0; l < max_kinematic_chain; l++) {            
            if (SegmentsKinematicChain(segment_indices(k, 0) - 1, l) == -1) break;
            
            //mexPrintf("l = %d\n", l);
            
            mwIndex joint_id = (mwIndex) SegmentsKinematicChain(segment_indices(k, 0) - 1, l) - 1;
            mwIndex segment_id = (mwIndex) JointsSegmentId(joint_id, 0) - 1;
            glm::vec4 u = glm::vec4(JointsAxis(joint_id, 0), JointsAxis(joint_id, 1), JointsAxis(joint_id, 2), 1);
            
            glm::dvec3 p = glm::dvec3(SegmentsGlobal(segment_id, 12), SegmentsGlobal(segment_id, 13), SegmentsGlobal(segment_id, 14));
            
            glm::mat4 T = glm::mat4(0);
            for (int x = 0; x < 4; x++) {
                for (int y = 0; y < 4; y++) {
                    T[y][x] = SegmentsGlobal(segment_id, 4 * y + x);
                    //mexPrintf("%f  ", T(i, j));
                }
                //mexPrintf("\n");
            }
            //mexPrintf("\n");
            
            // shape
            glm::vec4 w4  = T * glm::vec4(0, 1, 0, 1);
            glm::dvec3 w = glm::dvec3(w4[0] / w4[3], w4[1] / w4[3], w4[2] / w4[3]);
            w = w - p;
            
            double c = 1;
            bool last_in_kinematic_chain = SegmentsKinematicChain(segment_indices(k, 0) - 1, l + 1) == -1;
            if (l == 2 || last_in_kinematic_chain)
                c = glm::length(m - p) / beta(l, 0);
            glm::dvec3 cw = c * w;
            j.col(segment_id) = Eigen::Vector3d(cw[0], cw[1], cw[2]);
            
            // pose
            glm::vec4 v4 = T * u;
            glm::dvec3 v = glm::dvec3(v4[0] / v4[3], v4[1] / v4[3], v4[2] / v4[3]);
            v = v - p;
            glm::dvec3 vxd = glm::cross(v, m - p);
            j.col(num_segments - 1 + joint_id) = Eigen::Vector3d(vxd[0], vxd[1], vxd[2]);
            
            
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
            glm::dvec3 j_col = glm::dvec3(j(0, var), j(1, var), j(2, var));
            J(k, var) = - glm::dot(n, j_col);
            //mexPrintf("%f\t", J(k, var));
        }
        //mexPrintf("\n");
        F(k, 0) = glm::dot(n, d - m);
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
    
    mlx_array<mlx_double> DataPoints(mlx_size{num_points, 3}, in[3]);
    mlx_array<mlx_double> ModelPoints(mlx_size{num_points, 3}, in[4]);
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

