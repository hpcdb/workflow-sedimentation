/* 
 * File:   xdmf.h
 * Author: camata
 *
 * Created on February 26, 2015, 11:07 AM
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

class XDMF_IO
{
public:
    XDMF_IO(const Mesh& mesh, std::string basename);
    virtual ~XDMF_IO();
    string* write_timestep(EquationSystems &es, double time);
    int GetFileID() { return _timestep; }
    void SetFileID(int id) {_timestep = id;}
private:
    void libMesh_to_xdmf(std::vector<double> &coords, std::vector<int> &conn);
    std::map<dof_id_type, dof_id_type> node_map;
    
    std::vector<int>    g2l;
    std::string        _basename;
    int                _timestep;
    ElemType           _elemtype;
    int                _n_local_nodes;
    int                _n_local_elem;
    const Mesh &       _mesh;
    int                processor_id;
    int                n_processors;
        
    
    //
    void GetG2LMapping();
    void GetLocalCoordinates(std::vector<double> &coords);
    void GetElementConnectivity(std::vector<int> &conn);
    void GetLocalNavierStokesSolution(EquationSystems& es, std::vector<double> &velocity, std::vector<double> &pressure);
    void GetLocalSedimentationSolution(EquationSystems& es, std::vector<double> &sed, std::vector<double> &volume);
    
    void write_spatial_collection(double);
    void write_temporal_collection();
    
};

}

#endif	/* XDMF_H */

