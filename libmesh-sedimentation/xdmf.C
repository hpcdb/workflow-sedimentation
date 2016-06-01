/*
 * File:   xdmf.cpp
 * Author: camata
 *
 * Created on February 26, 2015, 11:07 AM
 */

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>
using namespace std;

#include "xdmf.h"

#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"


#include <hdf5.h>

void H5_WriteInteger(hid_t file_id, const char *datasetname, int *dataset, long size, bool using_compression)
{
    hid_t dset_id;
    hid_t dspace_id;
    hid_t dcpl_id;
    herr_t  status;
    hsize_t dims[1];
    dims[0] = size;

    dspace_id = H5Screate_simple(1,dims,NULL);
    if(using_compression)
    {
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        status  = H5Pset_chunk(dcpl_id,1,dims);
        status  = H5Pset_deflate(dcpl_id,6);
        dset_id = H5Dcreate2(file_id, datasetname, H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    } else {
        dset_id = H5Dcreate2(file_id,datasetname, H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset);
    if(using_compression) status = H5Pclose(dcpl_id);
    status = H5Dclose(dset_id);
    status = H5Sclose(dspace_id);

}

void H5_WriteDouble(hid_t file_id, const char *datasetname, double *dataset, long size, bool using_compression)
{
    hid_t dset_id;
    hid_t dspace_id;
    hid_t dcpl_id;
    herr_t  status;
    hsize_t dims[1];
    dims[0] = size;

    dspace_id = H5Screate_simple(1,dims,NULL);
    if(using_compression)
    {
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        status  = H5Pset_chunk(dcpl_id,1,dims);
        status  = H5Pset_deflate(dcpl_id,6);
        dset_id = H5Dcreate2(file_id, datasetname, H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    } else {
        dset_id = H5Dcreate2(file_id,datasetname, H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset);
    if(using_compression) status = H5Pclose(dcpl_id);
    status = H5Dclose(dset_id);
    status = H5Sclose(dspace_id);

}


namespace libMesh
{

XDMF_IO::XDMF_IO(const Mesh & mesh, std::string basename) : _mesh(mesh), _timestep(0),
                                                  _basename(basename)
{
    this->processor_id = libMesh::global_processor_id();
    this->n_processors = libMesh::global_n_processors();
}


void XDMF_IO::write_timestep(EquationSystems& es, double time)
{
    char filename[256];
    hid_t   file_id;
    herr_t  status;

    std::vector<double> velocity;
    std::vector<double> pressure;
    std::vector<double> coords;
    std::vector<int>    conn;
    std::vector<double> sediment;
    std::vector<double> volume;
    //printf("Criando arquivo hdf5: %s \n", filename);

    //this->GetG2LMapping();
    this->libMesh_to_xdmf(coords,conn);

    //this->GetLocalCoordinates(coords);
    //this->GetElementConnectivity(conn);

    this->GetLocalNavierStokesSolution(es,velocity, pressure);
    this->GetLocalSedimentationSolution(es,sediment,volume);

    //cout << "coords size = " << coords.size() << endl;

    sprintf(filename,"%s_%d_%03d_%05d.h5",this->_basename.c_str(),n_processors,processor_id, this->_timestep);



    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5_WriteInteger(file_id,"/conn",&conn[0],conn.size(),false);
    H5_WriteDouble(file_id,"/coords",&coords[0],coords.size(),false);
    H5_WriteDouble(file_id,"/velocity",&velocity[0],velocity.size(),false);
    H5_WriteDouble(file_id,"/pressure",&pressure[0],pressure.size(),false);
    H5_WriteDouble(file_id,"/sediment",&sediment[0],sediment.size(),false);
    H5_WriteDouble(file_id,"/volume" ,&volume[0],volume.size(),false);

    /* Close the file. */
    status = H5Fclose(file_id);

    write_spatial_collection(time);
    write_temporal_collection();

    this->_timestep++;
}


void XDMF_IO::write_spatial_collection(double time)
{


    if(processor_id == 0)
    {
        char filename[255];
        string elem_name;
        int nnoel = 0;
        if(this->_elemtype == HEX8) {
             elem_name= "Hexahedron";
             nnoel = 8;
        }
        else if(this->_elemtype == TET4) {
            elem_name = "Tetrahedron";
            nnoel = 4;
        }


        sprintf(filename,"%s_%d_%05d.xmf", this->_basename.c_str(), n_processors,this->_timestep);
        FILE * fxml = fopen(filename, "w");
        fprintf(fxml,"<?xml version=\"1.0\" ?>\n");
        fprintf(fxml,"<!-- DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [] --> \n");
        fprintf(fxml,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\"> \n");
        fprintf(fxml,"<Domain Name=\"libMesh\">\n");
        fprintf(fxml," <Grid Name=\"%s_%d_%05d\" GridType=\"Collection\" CollectionType=\"Spatial\">\n", this->_basename.c_str(), n_processors,this->_timestep );
        fprintf(fxml,"  <Time Type=\"Single\" Value=\"%f\" />\n", time);
        for(int p=0; p < n_processors; p++)
        {
            fprintf(fxml," <Grid Name=\"%s_%d_%03d\" Type=\"Uniform\"> \n",this->_basename.c_str(), n_processors, p);
            fprintf(fxml,"  <Topology Type=\"%s\" NumberOfElements=\"%d\">\n",elem_name.c_str(),this->_n_local_elem);
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\"> %s_%d_%03d_%05d.h5:/conn</DataItem>\n",this->_n_local_elem*nnoel,this->_basename.c_str(),n_processors,p,this->_timestep);
            fprintf(fxml,"  </Topology>\n");
            fprintf(fxml,"  <Geometry Type=\"XYZ\">\n");
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s_%d_%03d_%05d.h5:/coords</DataItem>\n",this->_n_local_nodes*3,this->_basename.c_str(),n_processors,p,this->_timestep);
            fprintf(fxml,"  </Geometry>\n");
            fprintf(fxml,"  <Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\">\n");
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s_%d_%03d_%05d.h5:/velocity</DataItem>\n",this->_n_local_nodes*3,this->_basename.c_str(),n_processors,p,this->_timestep);
            fprintf(fxml,"  </Attribute>\n");
            fprintf(fxml,"  <Attribute Name=\"pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s_%d_%03d_%05d.h5:/pressure</DataItem>\n",this->_n_local_nodes,this->_basename.c_str(),n_processors,p,this->_timestep);
            fprintf(fxml,"  </Attribute>\n");
            fprintf(fxml,"  <Attribute Name=\"sediment\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s_%d_%03d_%05d.h5:/sediment</DataItem>\n",this->_n_local_nodes,this->_basename.c_str(),n_processors,p,this->_timestep);
            fprintf(fxml,"  </Attribute>\n");
            fprintf(fxml,"  <Attribute Name=\"volume\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s_%d_%03d_%05d.h5:/volume</DataItem>\n",this->_n_local_nodes,this->_basename.c_str(),n_processors,p,this->_timestep);
            fprintf(fxml,"  </Attribute>\n");
            fprintf(fxml,"  </Grid>\n");
        }
        fprintf(fxml," </Grid>\n");
        fprintf(fxml,"</Domain>\n");
        fprintf(fxml,"</Xdmf>\n");
        fclose(fxml);
    }
}

void XDMF_IO::write_temporal_collection()
{
    if(processor_id == 0)
    {
        char filename[255];
        sprintf(filename,"%s_%d.xmf", this->_basename.c_str(), n_processors);
        FILE * fxml = fopen(filename, "w");
        fprintf(fxml,"<?xml version=\"1.0\" ?>\n");
        fprintf(fxml,"<!-- DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [] --> \n");
        fprintf(fxml,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\"> \n");
        fprintf(fxml,"<Domain Name=\"libMesh\">\n");
        fprintf(fxml," <Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
        for(int t=0; t <= this->_timestep; t++)
            fprintf(fxml,"  <xi:include href=\"%s_%d_%05d.xmf\" xpointer=\"xpointer(//Xdmf/Domain/Grid)\" />\n", this->_basename.c_str(), n_processors,t);
         fprintf(fxml," </Grid>\n");
        fprintf(fxml,"</Domain>\n");
        fprintf(fxml,"</Xdmf>\n");

        fclose(fxml);
    }
}

void XDMF_IO::libMesh_to_xdmf(std::vector<double>& coords, std::vector<int> &conn)
{
    node_map.clear();
    unsigned int local_node_counter = 0;
    MeshBase::const_node_iterator nd = _mesh.local_nodes_begin();
    MeshBase::const_node_iterator nd_end = _mesh.local_nodes_end();
    for (; nd != nd_end; nd++, ++local_node_counter)
    {
        Node* node = (*nd);
         // Fill mapping between global and local node numbers
        node_map[node->id()] = local_node_counter;
        double x = (*node)(0);
        double y = (*node)(1);
        double z = (*node)(2);
        coords.push_back(x);
        coords.push_back(y);
        coords.push_back(z);
    }

    MeshBase::const_element_iterator       it      = _mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end = _mesh.active_local_elements_end();

    unsigned active_element_counter = 0;
    for (; it != end; ++it, ++active_element_counter)
    {
        Elem *elem = *it;
        libmesh_assert(elem->type() == TET4 || elem->type() == HEX8);

        if(active_element_counter == 0 ) this->_elemtype = elem->type();
         for (unsigned int i=0; i<elem->n_nodes(); ++i)
         {
              dof_id_type global_node_id = elem->node(i);
              if (node_map.find(global_node_id) == node_map.end() ) {

                   const Node* node = _mesh.node_ptr(global_node_id);

                   // Error checking...
                    if (node == NULL)
                     libmesh_error_msg("Error getting pointer to node " << global_node_id << "!");

                    double x = (*node)(0);
                    double y = (*node)(1);
                    double z = (*node)(2);
                    coords.push_back(x);
                    coords.push_back(y);
                    coords.push_back(z);

                    node_map[global_node_id] = local_node_counter;

                    local_node_counter++;
              }

              conn.push_back(node_map[elem->node(i)]);
         }


    }

    this->_n_local_elem  = active_element_counter;
    this->_n_local_nodes = local_node_counter;

}


void XDMF_IO::GetG2LMapping()
{


    this->g2l.resize(_mesh.n_nodes());

    std::fill(this->g2l.begin(), this->g2l.end(), -1);

    MeshBase::const_element_iterator       it      = _mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_map = _mesh.active_local_elements_end();

    int n_nodes = 0;
    int n_elem = 0;
    for ( ; it != end_map; ++it)
    {
        Elem* elem = *it;
        if(n_nodes == 0 ) this->_elemtype = elem->type();
        for(unsigned int i =0; i < elem->n_nodes(); i++)
        {
            int id = elem->node(i);
            if(g2l[id] == -1) g2l[id] = n_nodes++;
        }

        libmesh_assert(elem->type() == this->_elemtype);
        n_elem++;
    }

    this->_n_local_nodes = n_nodes;
    this->_n_local_elem  = n_elem;

    //printf("n_local_nodes = %d\n", this->_n_local_nodes);
    //printf("n_local_elem  = %d\n", this->_n_local_elem );

}

void XDMF_IO::GetLocalCoordinates(std::vector<double>& coords)
{


    coords.resize(this->_n_local_nodes*3);

    MeshBase::const_element_iterator       it      = _mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_map = _mesh.active_local_elements_end();

    for ( ; it != end_map; ++it)
    {
        Elem* elem = *it;
        for(unsigned int i =0; i < elem->n_nodes(); i++)
        {
            int g_id           = elem->node(i);
            int l_id           = g2l[g_id];
            coords[l_id*3]     = elem->point(i)(0);
            coords[l_id*3+1]   = elem->point(i)(1);
            coords[l_id*3+2]   = elem->point(i)(2);
        }
    }

     //printf("coords \n");

}

void XDMF_IO::GetElementConnectivity(std::vector<int>& conn)
{


    MeshBase::const_element_iterator       it      = _mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_map = _mesh.active_local_elements_end();

    int n_elem_nodes;
    switch(this->_elemtype)
    {
        case TET4:
            n_elem_nodes = 4;
            break;
        case HEX8:
            n_elem_nodes = 8;
            break;
        default:
            n_elem_nodes = 8;
    }

    conn.resize(this->_n_local_elem*n_elem_nodes);

    int iel = 0;
    for ( ; it != end_map; ++it)
    {
        Elem* elem = *it;
        for(unsigned int i =0; i < elem->n_nodes(); i++)
        {
            int g_node_id = elem->node(i);
            int l_node_id = g2l[g_node_id];
            conn[iel*n_elem_nodes+i] = l_node_id;
        }
        iel++;
    }
     //printf("CONN \n");
}

void XDMF_IO::GetLocalNavierStokesSolution(EquationSystems& es, std::vector<double>& velocity, std::vector<double> &pressure)
{

  velocity.resize(this->_n_local_nodes*3);
  pressure.resize(this->_n_local_nodes);

    // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & navier_stokes_system =
    es.get_system<TransientLinearImplicitSystem> ("flow");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = navier_stokes_system.variable_number ("u");
  const unsigned int v_var = navier_stokes_system.variable_number ("v");
  const unsigned int w_var = navier_stokes_system.variable_number ("w");
  const unsigned int p_var = navier_stokes_system.variable_number ("p");
  DofMap & dof_map   = navier_stokes_system.get_dof_map();

  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_w, w_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);


      for(int i=0; i < elem->n_nodes(); i++)
      {
          int local_id = node_map[elem->node(i)];//g2l[elem->node(i)];
          velocity[local_id*3]   = navier_stokes_system.current_solution(dof_indices_u[i]);
          velocity[local_id*3+1] = navier_stokes_system.current_solution(dof_indices_v[i]);
          velocity[local_id*3+2] = navier_stokes_system.current_solution(dof_indices_w[i]);
          pressure[local_id]     = navier_stokes_system.current_solution(dof_indices_p[i]);
      }
    }

       //printf("Navier \n");

       //libMesh::CommWorld.barrier();


}

 void XDMF_IO::GetLocalSedimentationSolution(EquationSystems& es, std::vector<double>& sediment, std::vector<double>& volume)
{

  sediment.resize(this->_n_local_nodes);
  volume.resize(this->_n_local_nodes);

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // Get a reference to the Stokes system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("sediment");

  ExplicitSystem & deposition_system = es.get_system<ExplicitSystem>("deposition");
  //NumericVector<Number> &volume_vec = system.get_vector("volume");

  // Numeric ids corresponding to each variable in the system
  const unsigned int s_var            =  system.variable_number ("s");
  const unsigned int d_var            =  deposition_system.variable_number ("d");

  const DofMap & dof_map         = system.get_dof_map();
  const DofMap & dof_map_dep     = deposition_system.get_dof_map();

  std::vector<dof_id_type> dof_indices_s;
  std::vector<dof_id_type> dof_indices_d;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;

      dof_map.dof_indices    (elem, dof_indices_s ,  s_var);
      dof_map_dep.dof_indices(elem, dof_indices_d ,  d_var);

      for(int i=0; i < elem->n_nodes(); i++)
      {
          //int local_id = g2l[elem->node(i)];
          int local_id = node_map[elem->node(i)];
          sediment[local_id]    = system.current_solution(dof_indices_s[i]);
          volume[local_id]      = deposition_system.current_solution(dof_indices_d[i]);
      }
    }
}

 XDMF_IO::~XDMF_IO()
 {
     this->g2l.clear();
 }

}
