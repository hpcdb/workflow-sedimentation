/*
 * File:   xdmf.h
 * Author: camata
 *
 * Created on June 3, 2016, 2:13 PM
 */

#ifndef XDMF_H
#define	XDMF_H

#include <string>

#include "libmesh/libmesh_common.h"
#include "libmesh/mesh.h"

using namespace std;

namespace libMesh
{

class EquationSystems;
class MeshBase;
class System;

class XDMFWriter
{
    public:
        XDMFWriter(const Mesh& mesh);
        void set_file_name(std::string filename);
        void set_file_id(int n_time_file) {n_timestep = n_time_file; } 
        int  get_file_id() { return n_timestep; }
        string* write_time_step(EquationSystems &es, double time);
        virtual ~XDMFWriter();

    private:
        void libMesh_to_xdmf(std::vector<double> &coords, std::vector<int> &conn);
        void get_variable_solution(EquationSystems& es, int sys, int ivar, std::vector<double> &solution);
        void write_spatial_collection(EquationSystems& es, double time);
        void write_temporal_collection();

        // mapping between global and local ids
        std::map<dof_id_type, dof_id_type> g2l;
        std::string        basename;
        int                n_timestep;
        ElemType           elemtype;
        int                n_local_nodes;
        int                n_local_elem;
        const Mesh &       mesh;
        int                processor_id;
        int                n_processors;
};

}


#endif	/* XDMF_H */
