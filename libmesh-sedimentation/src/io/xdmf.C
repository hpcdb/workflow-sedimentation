
/* 
 * File:   xdmf.cpp
 * Author: camata
 * 
 * Created on February 26, 2015, 11:07 AM
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>
#include <unistd.h>
using namespace std;


#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"
#include "libmesh/fe_interface.h"


#ifdef LIBMESH_HAVE_HDF5
#include <hdf5.h>
#endif

#include "xdmf.h"

#define XDMF_FILE_SIZE 512

#ifdef LIBMESH_HAVE_HDF5
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
#else
void BIN_WriteInt(const char *datasetname, void *dataset, long size)
{
    FILE *fbin = fopen(datasetname,"w");
    fwrite(dataset,sizeof(int),size,fbin);
    fclose(fbin);
     
}

void BIN_WriteDouble(const char *datasetname, void *dataset, long size)
{
    FILE *fbin = fopen(datasetname,"w");
    fwrite(dataset,sizeof(double),size,fbin);
    fclose(fbin);
     
}

#endif


namespace libMesh
{

//
// XDMF Constructor. 
XDMFWriter::XDMFWriter(const Mesh & mesh) : mesh(mesh), 
   						n_timestep(0),
                                                basename("output"),
                                                dir("")
{
    this->processor_id = libMesh::global_processor_id();
    this->n_processors = libMesh::global_n_processors();
    this->is_mesh_changed = true;
    
   
}


XDMFWriter::~XDMFWriter()
{
    this->g2l.clear();
}

void XDMFWriter::set_file_name(std::string filename)
{
    this->basename = filename;
}

void XDMFWriter::set_dir_path(std::string path)
{
    this->dir = path;
    struct stat st = {0};
    if (stat(dir.c_str(), &st) == -1) {
        mkdir(dir.c_str(), 0700);
    }
}

string* XDMFWriter::write_time_step(EquationSystems& es, double time)
{
    char filename[XDMF_FILE_SIZE];
    char xdmf_filename[XDMF_FILE_SIZE];
    string* files = new string[2];
        
    std::vector<double> coords;
    std::vector<int>    conn;
    std::vector<double> solution;

    // Get local data
    this->libMesh_to_xdmf(coords,conn);
    
//#define DEBUG
#ifdef DEBUG
    int nnoel = 0;
     if(this->elemtype == HEX8) {
             nnoel = 8;
      }
        else if(this->elemtype == TET4 || this->elemtype == QUAD4) {
            
            nnoel = 4;
        } else if(this->elemtype == TRI3) {
             nnoel = 3;
        } 
    sprintf(filename,"%s_%d_%03d_%05d.dat",this->basename.c_str(),n_processors,processor_id, this->n_timestep);
    FILE * fout = fopen(filename, "w");
    fprintf(fout, "$NNODES %d\n", this->n_local_nodes);
    for(int i = 0; i < this->n_local_nodes; i++)
        fprintf(fout, "%8.8f %8.8f %8.8f\n", coords[i*3], coords[i*3+1], coords[i*3+2]);
    fprintf(fout, "$NELEM %d\n", this->n_local_elem);
   
    for(int i = 0; i < this->n_local_elem; i++) {
        for(int j = 0; j < nnoel; j++)
            fprintf(fout, "%d ", conn[i*nnoel + j]);
        fprintf(fout,"\n");  
    }
        
    fclose(fout);
    
#endif
    
    
    solution.resize(this->n_local_nodes);
    
    char stepdir[XDMF_FILE_SIZE];
    sprintf(stepdir,"%s/step%d", this->dir.c_str(), this->n_timestep);
    struct stat st = {0};
    if (stat(stepdir, &st) == -1) {
        mkdir(stepdir, 0700);
    }
    
    
#ifdef LIBMESH_HAVE_HDF5
    hid_t   file_id;
    herr_t  status;
    
    sprintf(filename,"%s/%s_%d_%03d_%05d.h5",stepdir,this->basename.c_str(),n_processors,processor_id, this->n_timestep);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5_WriteInteger(file_id,"/conn",&conn[0],conn.size(),false);
    H5_WriteDouble(file_id,"/coords",&coords[0],coords.size(),false);
    //LOOP over all system
    for(int ns=0; ns < es.n_systems(); ++ns)
    {
        
        for(int nv=0; nv < es.get_system(ns).n_vars(); ++nv )
        {
            this->get_variable_solution(es,ns,nv,solution);
            std::string dataset_name = "/"+es.get_system(ns).variable_name(nv);
            
            H5_WriteDouble(file_id,dataset_name.c_str(),&solution[0],solution.size(),false);
        }
    }
    
    status = H5Fclose(file_id);
#else

    sprintf(filename,"%s/%s.xyz.%d.%03d.%05d.bin",stepdir,this->basename.c_str(),n_processors,processor_id, this->n_timestep);
    BIN_WriteDouble(filename,&coords[0],coords.size());
    sprintf(filename,"%s/%s.con.%d.%03d.%05d.bin",stepdir,this->basename.c_str(),n_processors,processor_id, this->n_timestep);
    BIN_WriteInt(filename,&conn[0],conn.size());
    
    //LOOP over all system
    for(int ns=0; ns < es.n_systems(); ++ns)
    {
        
        for(int nv=0; nv < es.get_system(ns).n_vars(); ++nv )
        {
            this->get_variable_solution(es,ns,nv,solution);
            std::string v = es.get_system(ns).variable_name(nv);

            sprintf(filename,"%s/%s.%s.%d.%03d.%05d.bin",stepdir,this->basename.c_str(),v.c_str(),
                                                      n_processors,processor_id, this->n_timestep);
            BIN_WriteDouble(filename,&solution[0],solution.size());
        }
    }
    

#endif
    
    write_spatial_collection(es, time);
    write_temporal_collection();

    sprintf(xdmf_filename,"%s/%s_%d_%05d.xmf", this->dir.c_str(),this->basename.c_str(), n_processors,this->n_timestep);
    // char the_path[2048];
    // getcwd(the_path, 2048);
    files[0] = filename;
    files[1] = xdmf_filename;
  
    this->n_timestep++;
    return files;
}


void XDMFWriter::write_spatial_collection(EquationSystems& es, double time)
{
    int data[2];
    int * local_data = new int [2*this->n_processors];
    
    data[0] = this->n_local_nodes;
    data[1] = this->n_local_elem;
    MPI_Gather(data, 2, MPI_INT, local_data, 2, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(processor_id == 0)
    {
        char filename[255];
        string elem_name;
        int nnoel = 0;
        if(this->elemtype == HEX8) {
             elem_name= "Hexahedron";
             nnoel = 8;
        }
        else if(this->elemtype == TET4) {
            elem_name = "Tetrahedron";
            nnoel = 4;
        } else if(this->elemtype == TRI3) {
             elem_name= "Triangle";
             nnoel = 3;
        }
        else if(this->elemtype == QUAD4) {
            elem_name = "Quadrilateral";
            nnoel = 4;
        } 
        
        char stepdir[XDMF_FILE_SIZE];
        sprintf(stepdir,"%s/step%d", this->dir.c_str(), this->n_timestep);
        char stepfile[XDMF_FILE_SIZE];
        sprintf(stepfile,"step%d/%s",this->n_timestep ,this->basename.c_str());
        
       
        sprintf(filename,"%s/%s_%d_%05d.xmf", this->dir.c_str(),this->basename.c_str(), n_processors,this->n_timestep);
        FILE * fxml = fopen(filename, "w");
        fprintf(fxml,"<?xml version=\"1.0\" ?>\n");
        fprintf(fxml,"<!-- DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [] --> \n");
        fprintf(fxml,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\"> \n");
        fprintf(fxml,"<Domain Name=\"libMesh\">\n");
        fprintf(fxml," <Grid Name=\"%s_%d_%05d\" GridType=\"Collection\" CollectionType=\"Spatial\">\n", stepfile, n_processors,this->n_timestep );
        fprintf(fxml,"  <Time Type=\"Single\" Value=\"%f\" />\n", time);
        for(int p=0; p < n_processors; p++)
        {
            int lnodes = local_data[p*2];
            int lelem  = local_data[p*2+1];
            fprintf(fxml," <Grid Name=\"%s_%d_%03d\" Type=\"Uniform\"> \n",this->basename.c_str(), n_processors, p);  
            fprintf(fxml,"  <Topology Type=\"%s\" NumberOfElements=\"%d\"  BaseOffset=\"0\">\n",elem_name.c_str(),lelem);
#ifdef LIBMESH_HAVE_HDF5
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\"> %s_%d_%03d_%05d.h5:/conn</DataItem>\n",lelem*nnoel,stepfile,n_processors,p,this->n_timestep);
#else
           fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Int\" Format=\"Binary\" Endian=\"Little\"> %s.con.%d.%03d.%05d.bin </DataItem>\n",lelem*nnoel,stepfile,n_processors,p,this->n_timestep);
#endif

            fprintf(fxml,"  </Topology>\n");
            fprintf(fxml,"  <Geometry Type=\"XYZ\">\n");
#ifdef LIBMESH_HAVE_HDF5
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">%s_%d_%03d_%05d.h5:/coords</DataItem>\n",lnodes*3,stepfile,n_processors,p,this->n_timestep);
#else
            fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">%s.xyz.%d.%03d.%05d.bin </DataItem>\n",lnodes*3,stepfile,n_processors,p,this->n_timestep);
#endif
            fprintf(fxml,"  </Geometry>\n");
            
            //LOOP over all system
            for(int ns=0; ns < es.n_systems(); ++ns)
            {
        
                for(int nv=0; nv < es.get_system(ns).n_vars(); ++nv )
                {
                    std::string dataset_name = es.get_system(ns).variable_name(nv);
                    fprintf(fxml,"  <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", dataset_name.c_str());
#ifdef LIBMESH_HAVE_HDF5
                    fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\"> %s_%d_%03d_%05d.h5:/%s</DataItem>\n",lnodes,stepfile,n_processors,p,this->n_timestep, dataset_name.c_str());
#else
		   fprintf(fxml,"    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"Binary\" Endian=\"Little\">%s.%s.%d.%03d.%05d.bin </DataItem>\n",lnodes,stepfile,dataset_name.c_str(),n_processors,p,this->n_timestep);
#endif
                    fprintf(fxml,"  </Attribute>\n");
            
                }
           }
        
           fprintf(fxml,"  </Grid>\n");
        }
        fprintf(fxml," </Grid>\n");
        fprintf(fxml,"</Domain>\n"); 
        fprintf(fxml,"</Xdmf>\n");
        fclose(fxml);
        
        delete [] local_data;
    }
}

void XDMFWriter::write_temporal_collection()
{
    if(processor_id == 0)
    {
        char filename[XDMF_FILE_SIZE];   
        sprintf(filename,"%s/%s_%d.xmf", this->dir.c_str(),this->basename.c_str(), n_processors);
        FILE * fxml = fopen(filename, "w");
        fprintf(fxml,"<?xml version=\"1.0\" ?>\n");
        fprintf(fxml,"<!-- DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" [] --> \n");
        fprintf(fxml,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\"> \n");
        fprintf(fxml,"<Domain Name=\"libMesh\">\n");
        fprintf(fxml," <Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
        for(int t=0; t <= this->n_timestep; t++)
            fprintf(fxml,"  <xi:include href=\"%s_%d_%05d.xmf\" xpointer=\"xpointer(//Xdmf/Domain/Grid)\" />\n", this->basename.c_str(), n_processors,t);
         fprintf(fxml," </Grid>\n");
        fprintf(fxml,"</Domain>\n"); 
        fprintf(fxml,"</Xdmf>\n");
        
        fclose(fxml);
    }
}


void XDMFWriter::mesh_changed_on()
{
    this->is_mesh_changed = true;
}


void XDMFWriter::libMesh_to_xdmf(std::vector<double>& coords, std::vector<int> &conn)
{
    
    //if(!this->is_mesh_changed) return;
    
    g2l.clear();
    coords.clear();
    conn.clear();
 
    // NEW 
    /*
    dof_id_type local_node_counter = 0;
    MeshBase::const_node_iterator n_it        = mesh.local_nodes_begin();
    const MeshBase::const_node_iterator n_end = mesh.local_nodes_end();

    for (; n_it != n_end; ++n_it)
    {
        const Node* node = *n_it;
        dof_id_type gid = node->id();
        g2l.insert(std::pair<int, int>(gid,local_node_counter));
        //g2l[gid] = local_node_counter;
        double x = (*node)(0);
        double y = (*node)(1);
        double z = (*node)(2);
        coords.push_back(x);
        coords.push_back(y);
        coords.push_back(z);
        local_node_counter++;
       
    }
    */
   
   
  
   //OLD  
    MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
   
    
    unsigned int local_node_counter = 0;
    for (; it != end; ++it)
    {
        const Elem *elem = *it;
        libmesh_assert(elem->type() == TET4 || elem->type() == HEX8 || 
                       elem->type() == TRI3 || elem->type() == QUAD4);
        
        for (unsigned int i=0; i<elem->n_nodes(); ++i)
        {
              dof_id_type global_node_id = elem->node(i);
              if (g2l.find(global_node_id) == g2l.end() ) 
              {
                  
                    double x = elem->point(i)(0);
                    double y = elem->point(i)(1);
                    double z = elem->point(i)(2);
                    coords.push_back(x);
                    coords.push_back(y);
                    coords.push_back(z);
                    
                    g2l[global_node_id] = local_node_counter;
                    
                    local_node_counter++;
              }
         }
        
    }
    

    libmesh_assert(local_node_counter == mesh.n_local_nodes());
    
    unsigned int local_elem_counter = 0;
    
    it  = mesh.active_local_elements_begin();
    for (; it != end; ++it)
     {
        const Elem *elem = *it;
        if(local_elem_counter == 0) this->elemtype = elem->type();
        for (unsigned int i=0; i<elem->n_nodes(); ++i)
        {
            
            libmesh_assert(g2l.find(elem->node(i)) != g2l.end());
            int local_node = g2l[elem->node(i)];
            libmesh_assert(local_node < local_node_counter);
            conn.push_back(local_node);
        }
        local_elem_counter++;
        
    }
     
    libmesh_assert(local_elem_counter == mesh.n_active_local_elem());
    
    this->n_local_elem  = local_elem_counter;
    this->n_local_nodes = local_node_counter;
    this->is_mesh_changed =false;
   
}
    

void XDMFWriter::get_variable_solution(EquationSystems& es, int sys, int ivar, std::vector<double> &solution)
{
    const unsigned int dim = mesh.mesh_dimension();
    const System & system  = es.get_system(sys);
    const DofMap & dof_map = system.get_dof_map();
 
    std::vector<dof_id_type> dof_indices;
    std::vector<double>       nodal_soln;
    std::vector<double>       elem_soln;
    
    const FEType & fe_type    = system.variable_type(ivar);
    
    NumericVector<Number> & sys_soln(*system.current_local_solution);
    //std::cout << "Solution size: " << sys_soln.size() << std::endl;
    //std::cout << "# local Nodes: " << mesh.n_local_nodes() << std::endl;
    
    /*
    MeshBase::const_node_iterator it        = mesh.local_nodes_begin();
    const MeshBase::const_node_iterator end = mesh.local_nodes_end();
    for (; it != end; ++it)
    {
        const Node* node = *it;
        
        dof_id_type local_id       = g2l[node->id()];
        dof_id_type dof_id         = node->dof_number(sys,ivar,0);
        solution[local_id] = sys_soln(dof_id);
        
    } **/
   
    
    MeshBase::const_element_iterator       it  = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
    
    for (; it != end; ++it)
    {
        const Elem *elem = *it;
        
        dof_map.dof_indices (elem, dof_indices, ivar);
        elem_soln.resize(dof_indices.size());
        for(int i = 0; i < dof_indices.size(); ++i)
            elem_soln[i] = sys_soln(dof_indices[i]);
        
        
        FEInterface::nodal_soln (dim,fe_type, elem, elem_soln,nodal_soln);
        
        libmesh_assert_equal_to (nodal_soln.size(), elem->n_nodes());
        for (unsigned int n=0; n<elem->n_nodes(); n++) {
            int local_id = g2l[elem->node(n)];
            solution[local_id] = nodal_soln[n];
        }
    }
    
  
}

}
