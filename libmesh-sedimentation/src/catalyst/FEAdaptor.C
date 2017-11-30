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
#include "vtkCPPythonAdaptorAPI.h"

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

std::pair<int, int> CellTypeLibMeshToVTK(ElemType etype) {
    switch (etype) {
        case TRI3:
            return std::pair<int, int> (VTK_TRIANGLE, 3);
        case QUAD4:
            return std::pair<int, int> (VTK_QUAD, 4);
        case TET4:
            return std::pair<int, int> (VTK_TETRA, 4);
        case HEX8:
            return std::pair<int, int> (VTK_HEXAHEDRON, 8);
        default:
            std::cout << "Element type not supported" << endl;
            return std::pair<int, int>(-1, -1);
    }

}


namespace {

    vtkCPProcessor *Processor    = NULL;
    vtkUnstructuredGrid *VTKGrid = NULL;
    int rebuild_grid             = false;

    void CreateVTKGrid(EquationSystems &eq, std::map<unsigned int, unsigned int> &g2l) {     
        const MeshBase &mesh = eq.get_mesh();
        g2l.clear();

        std::cout << "  Creating Grid..." << std::endl;

        int unsigned ncells  = mesh.n_local_elem();
        int unsigned npoints = mesh.n_local_nodes();

        // create the points information
        vtkNew<vtkDoubleArray> pointArray;
        pointArray->SetNumberOfComponents(3);
        //pointArray->Allocate(static_cast<vtkIdType> (npoints * 3));


        MeshBase::const_element_iterator e_iter = mesh.active_local_elements_begin();
        MeshBase::const_element_iterator e_end  = mesh.active_local_elements_end();
        ElemType etype = (*e_iter)->type();
        int nodal_counter = 0;
        for (; e_iter != e_end; e_iter++) {
            const Elem* elem = *e_iter;
            for (int n = 0; n < elem->n_nodes(); n++) {
                int g_id = elem->node(n);
                if (g2l.find(g_id) == g2l.end()) {
                    g2l[g_id] = nodal_counter;

                    pointArray->InsertNextTuple3(elem->point(n)(0),
                            elem->point(n)(1),
                            elem->point(n)(2));
                    nodal_counter++;
                }
            }
        }

        vtkNew<vtkPoints> points;
        points->SetData(pointArray.GetPointer());
        VTKGrid->SetPoints(points.GetPointer());

        libmesh_assert(g2l.size() == npoints);


        e_iter = mesh.active_local_elements_begin();
        for (; e_iter != e_end; ++e_iter) {
            vtkIdType tmp[8];
            const Elem* elem = *e_iter;
            std::pair<int, int> elemmap = CellTypeLibMeshToVTK(elem->type());
            for (int ino = 0; ino < elem->n_nodes(); ++ino) {
                tmp[ino] = g2l[elem->node(ino)];
            }

            if (elemmap.first == VTK_TETRA) {
                int t = tmp[0];
                tmp[0] = tmp[1];
                tmp[1] = t;
            }

            VTKGrid->InsertNextCell(elemmap.first, elemmap.second, tmp);

        }

    }

    void get_variable_solution(const System& system, int var_number, std::vector<double> &solution, std::map<unsigned int, unsigned int> &g2l) {

        const MeshBase &mesh = system.get_mesh();
        const unsigned int dim = mesh.mesh_dimension();
        const DofMap & dof_map = system.get_dof_map();

        std::vector<dof_id_type> dof_indices;
        std::vector<double> nodal_soln;
        std::vector<double> elem_soln;

        const FEType & fe_type = system.variable_type(var_number);

        NumericVector<Number> & sys_soln(*system.current_local_solution);

        MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator end = mesh.active_local_elements_end();

        for (; it != end; ++it) {

            const Elem *elem = *it;

            dof_map.dof_indices(elem, dof_indices, var_number);

            elem_soln.resize(dof_indices.size());
            for (int i = 0; i < dof_indices.size(); ++i)
                elem_soln[i] = sys_soln(dof_indices[i]);

            FEInterface::nodal_soln(dim, fe_type, elem, elem_soln, nodal_soln);

            libmesh_assert_equal_to(nodal_soln.size(), elem->n_nodes());

            for (unsigned int n = 0; n < elem->n_nodes(); n++) {
                int local_id = g2l[elem->node(n)];
                solution[local_id] = nodal_soln[n];
            }
        }

    }

    void UpdateFields(EquationSystems &eq,
            std::map<unsigned int, unsigned int> & libmesh_global_to_local_map,
            vtkCPDataDescription* dataDescription) {
        vtkCPInputDataDescription* idd = dataDescription->GetInputDescriptionByName("input");

        if (VTKGrid == NULL) {
            VTKGrid = vtkUnstructuredGrid::New();
            CreateVTKGrid(eq, libmesh_global_to_local_map);
        }

        // If AMR is used, we need to rebuild the VTKGrid.
        if (VTKGrid != NULL && rebuild_grid) {
            libmesh_global_to_local_map.clear();
            VTKGrid->Delete();
            VTKGrid = vtkUnstructuredGrid::New();
            CreateVTKGrid(eq, libmesh_global_to_local_map);
            rebuild_grid = false;

        }

        std::cout << "  Updating Fields Data..." << std::endl;

        const MeshBase & mesh = eq.get_mesh();

        vtkIdType NumberOfNodes = static_cast<vtkIdType> (VTKGrid->GetNumberOfPoints());


        unsigned int n_systems = eq.n_systems();

        // now add numerical fields data
        if (VTKGrid->GetPointData()->GetNumberOfArrays() == 0) {

            for (int s = 0; s < n_systems; s++) {
                // Getting the system
                const System & sys = eq.get_system(s);

                for (int v = 0; v < sys.n_vars(); v++) {
                    const std::string var_name = sys.variable_name(v);

                    if (idd->IsFieldNeeded(var_name.c_str())) {

                        vtkDoubleArray* pointData = vtkDoubleArray::New();
                        pointData->SetName(var_name.c_str());
                        pointData->SetNumberOfComponents(1);
                        pointData->SetNumberOfTuples(NumberOfNodes);
                        VTKGrid->GetPointData()->AddArray(pointData);
                        pointData->Delete();

                    }
                }

            }

        }

        {
            std::vector<double> solution;
            solution.resize(libmesh_global_to_local_map.size());


            for (int s = 0; s < n_systems; s++) {
                const System & system = eq.get_system(s);

                for (int v = 0; v < system.n_vars(); v++) {
                    const std::string var_name = system.variable_name(v);
                    if (idd->IsFieldNeeded(var_name.c_str())) {

                        get_variable_solution(system, v, solution, libmesh_global_to_local_map);
                        vtkDoubleArray* data = vtkDoubleArray::SafeDownCast(VTKGrid->GetPointData()->GetArray(var_name.c_str()));
                        for (int i = 0; i < VTKGrid->GetNumberOfPoints(); i++)
                            data->SetTuple(i, &solution[i]);

                    }
                }

            }
        }
    }

}

namespace FEAdaptor {

    void mark_to_rebuild_grid() {
        rebuild_grid = true;
    }

    std::map<unsigned int, unsigned int> libmesh_global_to_local_map;

    void Initialize(int numScripts, string extractionScript, vector<string> visualizationScripts) {

        if (Processor == NULL) {
            Processor = vtkCPProcessor::New();
            Processor->Initialize();
        } else {
            Processor->RemoveAllPipelines();
        }
        
        if (numScripts > 0 && (extractionScript != "not defined")) {
            vtkNew<vtkCPPythonScriptPipeline> extraction;
            extraction->Initialize(extractionScript.c_str());
            Processor->AddPipeline(extraction.GetPointer());
        }

        for (int i = 0; i < visualizationScripts.size(); i++) {
            if (numScripts > 0 && (visualizationScripts[i] != "not defined")) {
                vtkNew<vtkCPPythonScriptPipeline> visualization;
                visualization->Initialize(visualizationScripts[i].c_str());
                Processor->AddPipeline(visualization.GetPointer());
            }
        }

    }

    void Finalize() {
        if (Processor) {
            Processor->Delete();
            Processor = NULL;
        }
        if (VTKGrid) {
            VTKGrid->Delete();
            VTKGrid = NULL;
        }
    }

    void CoProcess(EquationSystems &eq, double time, unsigned int timeStep, unsigned int analysisInterval, bool lastTimeStep = false, bool using_amr = false) {
        vtkNew<vtkCPDataDescription> dataDescription;
        dataDescription->AddInput("input");
        dataDescription->SetTimeData(time, timeStep);

        if (lastTimeStep == true) {
            // assume that we want to all the pipelines to execute if it
            // is the last time step.
            dataDescription->ForceOutputOn();
        }
        if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {
            UpdateFields(eq, libmesh_global_to_local_map, dataDescription.GetPointer());
            dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
            dataDescription->ForceOutputOn();
            Processor->CoProcess(dataDescription.GetPointer());
        }
    }

} // end of Catalyst namespace

#endif