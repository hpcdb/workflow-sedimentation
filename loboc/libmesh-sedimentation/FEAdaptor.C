#include <iostream>
#include <map>
#include "FEAdaptor.h"

#ifdef USE_CATALYST

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCPAdaptorAPI.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkMultiProcessController.h>

#include "libmesh/libmesh.h"


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

#include <iostream>
#include <string>

using namespace std;

namespace
{

  vtkCPProcessor* Processor    = NULL;
  vtkUnstructuredGrid* VTKGrid = NULL;
  
  void BuildVTKGrid(EquationSystems &eq, std::map<int, int> &g2l)
  {
    
    //std::vector<double> coords;
    g2l.clear();
    
    //std::cout << "BuildVTKGrid BEGIN" << std::endl;
    const MeshBase &mesh = eq.get_mesh();
    
    int numberOfCells  = mesh.n_local_elem();
    int numberOfPoints = mesh.n_local_nodes(); 
    
    vtkNew<vtkDoubleArray> pointArray;
   
    pointArray->SetNumberOfComponents(3);
    pointArray->Allocate(static_cast<vtkIdType>(numberOfPoints*3));

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    ElemType etype = (*el)->type();
    int counter = 0;
    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;


      for(int ino = 0; ino <  elem->n_nodes(); ++ino)
      {
          int g_id = elem->node(ino);
          
          if(g2l.find(g_id) == g2l.end() )
          {
              double p[3];
              g2l[g_id] = counter;
              p[0] = elem->point(ino)(0);
              p[1] = elem->point(ino)(1);
              p[2] = elem->point(ino)(2);
              pointArray->InsertNextTuple(p);
              counter++;
          }
      }
    }
  
    vtkNew<vtkPoints> points;
    points->SetData(pointArray.GetPointer());
    VTKGrid->SetPoints(points.GetPointer());
    
    //std::cout << "VTK INSERT NODES" << std::endl;
    
    libmesh_assert(g2l.size() == numberOfPoints);
  
   
    // create the cells
     int np = 0;
     int vtk_type;
    if(etype == HEX8) {
       VTKGrid->Allocate(static_cast<vtkIdType>(numberOfCells*9));
       np = 8;
       vtk_type = VTK_HEXAHEDRON;

     }
     else if (etype == TET4) {
       VTKGrid->Allocate(static_cast<vtkIdType>(numberOfCells*5));
        np =4;
        vtk_type = VTK_TETRA;
     }
    
    el     = mesh.active_local_elements_begin();
    for ( ; el != end_el; ++el)
    {
        vtkIdType tmp[8];
        const Elem* elem = *el;
        for(int ino = 0; ino < elem->n_nodes(); ++ino)
        {
            tmp[ino] = g2l[elem->node(ino)];
        }

         VTKGrid->InsertNextCell(vtk_type,np,  tmp);
               
    }
    
    //std::cout << "VTK INSERT ELEM" << std::endl;
    
   
  }
  
  void get_variable_solution(EquationSystems& es, int sys, int ivar, std::vector<double> &solution, std::map<int, int> &g2l)
{
    const MeshBase &mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const System & system  = es.get_system(sys);
    const DofMap & dof_map = system.get_dof_map();
 
    std::vector<dof_id_type> dof_indices;
    std::vector<double>       nodal_soln;
    std::vector<double>       elem_soln;
    
    const FEType & fe_type    = system.variable_type(ivar);
    
    NumericVector<Number> & sys_soln(*system.current_local_solution);
    
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
  
  void UpdateVTKAttributes(EquationSystems &es, std::map<int, int> &g2l )
  {
     //std::cout << "UpdateVTKAttributes BEGIN" << std::endl;
     std::vector<double> solution;
     solution.resize(g2l.size());
     

     if(VTKGrid->GetPointData()->GetNumberOfArrays() == 0)
      {
         for(int ns=0; ns < es.n_systems(); ++ns) {
            for(int nv=0; nv < es.get_system(ns).n_vars(); ++nv ) 
            {
                 
                get_variable_solution(es,ns,nv,solution,g2l);
                vtkNew<vtkDoubleArray> data;
                data->SetName(es.get_system(ns).variable_name(nv).c_str());
                data->SetNumberOfComponents(1);
                data->SetNumberOfTuples(static_cast<vtkIdType>(VTKGrid->GetNumberOfPoints()));
                for(int i = 0; i < VTKGrid->GetNumberOfPoints(); i++)
                     data->SetTuple(i,&solution[i]);
                 VTKGrid->GetPointData()->AddArray(data.GetPointer());
            }
         }
     } 
     else
     {
       for(int ns=0; ns < es.n_systems(); ++ns) {
        for(int nv=0; nv < es.get_system(ns).n_vars(); ++nv )
        {
            get_variable_solution(es,ns,nv,solution,g2l);
            std::string var_name = es.get_system(ns).variable_name(nv);
            vtkDoubleArray* data = vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray(var_name.c_str()));
            for(int i = 0; i < VTKGrid->GetNumberOfPoints(); i++)
                     data->SetTuple(i,&solution[i]);
       
        }
      }
    }
     
  }
  
  void BuildVTKDataStructures(EquationSystems &eq, std::map<int, int> &g2l, bool using_amr)
  {
      
    
    //std::cout << "BuildVTKDataStructures BEGIN" << std::endl;
    if(VTKGrid == NULL)
    {
        // The grid structure isn't changing so we only build it
        // the first time it's needed. If we needed the memory
        // we could delete it and rebuild as necessary.
        VTKGrid = vtkUnstructuredGrid::New();
        BuildVTKGrid(eq, g2l);
    }
    
    if(VTKGrid != NULL && using_amr)
    {
        g2l.clear();
        VTKGrid->Delete();
        VTKGrid = vtkUnstructuredGrid::New();
        BuildVTKGrid(eq, g2l);
        
    }
    
    UpdateVTKAttributes(eq, g2l);
    
  }

}

namespace FEAdaptor
{

  std::map<int, int>                   g2l;
  int numScripts;
  vtkNew<vtkCPPythonScriptPipeline> extraction;
  vtkNew<vtkCPPythonScriptPipeline> visualization;

  
  void Initialize(int numScripts, string extractionScript, string visualizationScript)
  {
      
    //std::cout << "COPROCESSING INIT BEGIN" << std::endl;
    if(Processor == NULL)
      {
        Processor = vtkCPProcessor::New();
        Processor->Initialize();
      }
    else
      {
        Processor->RemoveAllPipelines();
      }

      if(numScripts > 1){
        extraction->Initialize(extractionScript.c_str());
        Processor->AddPipeline(extraction.GetPointer());
      }

      if(numScripts > 2){
        visualization->Initialize(visualizationScript.c_str());
        Processor->AddPipeline(visualization.GetPointer());
      }
     //std::cout << "COPROCESSING INIT END" << std::endl;
  }

  void Finalize()
  {
    if(Processor)
      {
        Processor->Delete();
        Processor = NULL;
      }
    if(VTKGrid)
      {
         VTKGrid->Delete();
         VTKGrid = NULL;
      }
  }

  void CoProcess(int numScripts, string extractionScript, string visualizationScript, EquationSystems &eq, double time, unsigned int timeStep, 
    bool lastTimeStep = false, bool using_amr = false)
  {
    //std::cout << "COPROCESSING BEGIN" << std::endl;
    if(numScripts > 1){
      Processor->RemovePipeline(extraction.GetPointer());
      extraction->Initialize(extractionScript.c_str());
      Processor->AddPipeline(extraction.GetPointer());
    }

    vtkNew<vtkCPDataDescription> dataDescription;
    dataDescription->AddInput("input");
    dataDescription->SetTimeData(time, timeStep);
    if(lastTimeStep == true)
    {
      // assume that we want to all the pipelines to execute if it
      // is the last time step.
      dataDescription->ForceOutputOn();
    }
    if(Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
    {
      BuildVTKDataStructures(eq, g2l, using_amr);
      dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
      dataDescription->ForceOutputOn();
      Processor->CoProcess(dataDescription.GetPointer());
    }
  }
  
} // end of Catalyst namespace

#endif