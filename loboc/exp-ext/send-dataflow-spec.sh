#!/bin/bash
SIMULATION_DIR=`pwd`
curl -X POST http://localhost:22000/pde/dataflow/json \
    -H "Content-Type: application/text" \
    --data '{
  "transformations": [
    {
      "sets": [
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            }
          ],
          "tag": "iinputmesh",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "dim",
              "type": "NUMERIC"
            },
            {
              "name": "mesh_file",
              "type": "TEXT"
            },
            {
              "name": "restart_control",
              "type": "TEXT"
            }
          ],
          "tag": "oinputmesh",
          "type": "OUTPUT"
        }
      ],
      "tag": "inputmesh",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::InputMesh"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "inputmesh",
          "tag": "oinputmesh",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "r_fraction",
              "type": "NUMERIC"
            },
            {
              "name": "c_fraction",
              "type": "NUMERIC"
            },
            {
              "name": "max_h_level",
              "type": "NUMERIC"
            },
            {
              "name": "hlevels",
              "type": "NUMERIC"
            },
            {
              "name": "first_step_refinement",
              "type": "TEXT"
            },
            {
              "name": "amrc_flow_transp",
              "type": "TEXT"
            },
            {
              "name": "ref_interval",
              "type": "NUMERIC"
            },
            {
              "name": "max_r_steps",
              "type": "NUMERIC"
            }
          ],
          "tag": "oamrconfig",
          "type": "OUTPUT"
        }
      ],
      "tag": "amrconfig",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::AMRConfig"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "inputmesh",
          "tag": "oinputmesh",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "reynolds",
              "type": "NUMERIC"
            },
            {
              "name": "gr",
              "type": "NUMERIC"
            },
            {
              "name": "sc",
              "type": "NUMERIC"
            },
            {
              "name": "us",
              "type": "NUMERIC"
            },
            {
              "name": "diffusivity",
              "type": "NUMERIC"
            },
            {
              "name": "xlock",
              "type": "NUMERIC"
            },
            {
              "name": "fopc",
              "type": "NUMERIC"
            },
            {
              "name": "theta",
              "type": "NUMERIC"
            },
            {
              "name": "ex",
              "type": "NUMERIC"
            },
            {
              "name": "ey",
              "type": "NUMERIC"
            },
            {
              "name": "ez",
              "type": "NUMERIC"
            },
            {
              "name": "c_factor",
              "type": "NUMERIC"
            }
          ],
          "tag": "ocreateequationsystems",
          "type": "OUTPUT"
        }
      ],
      "tag": "createequationsystems",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::CreateEquationSystems"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "inputmesh",
          "tag": "oinputmesh",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "model_name",
              "type": "TEXT"
            },
            {
              "name": "dt_min",
              "type": "NUMERIC"
            },
            {
              "name": "dt_max",
              "type": "NUMERIC"
            },
            {
              "name": "tol_u",
              "type": "NUMERIC"
            },
            {
              "name": "tol_s",
              "type": "NUMERIC"
            },
            {
              "name": "kp",
              "type": "NUMERIC"
            },
            {
              "name": "ki",
              "type": "NUMERIC"
            },
            {
              "name": "kd",
              "type": "NUMERIC"
            },
            {
              "name": "nsa_max",
              "type": "NUMERIC"
            },
            {
              "name": "nsa_target_flow",
              "type": "NUMERIC"
            },
            {
              "name": "nsa_target_transport",
              "type": "NUMERIC"
            },
            {
              "name": "nsa_limit_flow",
              "type": "NUMERIC"
            },
            {
              "name": "nsa_limit_transport",
              "type": "NUMERIC"
            },
            {
              "name": "mult_factor_max",
              "type": "NUMERIC"
            },
            {
              "name": "mult_factor_min",
              "type": "NUMERIC"
            },
            {
              "name": "pc11_theta",
              "type": "NUMERIC"
            },
            {
              "name": "alpha",
              "type": "NUMERIC"
            },
            {
              "name": "k_exp",
              "type": "NUMERIC"
            },
            {
              "name": "s_min",
              "type": "NUMERIC"
            },
            {
              "name": "s_max",
              "type": "NUMERIC"
            },
            {
              "name": "reduct_factor",
              "type": "NUMERIC"
            },
            {
              "name": "complete_flow_norm",
              "type": "TEXT"
            }
          ],
          "tag": "otimestepcontrolconfig",
          "type": "OUTPUT"
        }
      ],
      "tag": "timestepcontrolconfig",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::TimeStepControlConfig"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "inputmesh",
          "tag": "oinputmesh",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "dpath",
              "type": "TEXT"
            },
            {
              "name": "rname",
              "type": "TEXT"
            },
            {
              "name": "write_interval",
              "type": "NUMERIC"
            },
            {
              "name": "catalyst_interval",
              "type": "NUMERIC"
            },
            {
              "name": "write_restart",
              "type": "TEXT"
            }
          ],
          "tag": "oioconfig",
          "type": "OUTPUT"
        }
      ],
      "tag": "ioconfig",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::IOConfig"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "createequationsystems",
          "tag": "ocreateequationsystems",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "dt",
              "type": "NUMERIC"
            },
            {
              "name": "tmax",
              "type": "NUMERIC"
            },
            {
              "name": "n_time_steps",
              "type": "NUMERIC"
            },
            {
              "name": "n_nonlinear_steps",
              "type": "NUMERIC"
            },
            {
              "name": "nonlinear_tolerance",
              "type": "NUMERIC"
            },
            {
              "name": "max_linear_iters",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "TEXT"
            }
          ],
          "tag": "ogetmaximumiterationstoflow",
          "type": "OUTPUT"
        }
      ],
      "tag": "getmaximumiterationstoflow",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::getMaximumIterationsToFlow"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstoflow",
          "tag": "ogetmaximumiterationstoflow",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "dt",
              "type": "NUMERIC"
            },
            {
              "name": "tmax",
              "type": "NUMERIC"
            },
            {
              "name": "n_time_steps",
              "type": "NUMERIC"
            },
            {
              "name": "n_nonlinear_steps",
              "type": "NUMERIC"
            },
            {
              "name": "nonlinear_tolerance",
              "type": "NUMERIC"
            },
            {
              "name": "max_linear_iters",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "TEXT"
            }
          ],
          "tag": "ogetmaximumiterationstotransport",
          "type": "OUTPUT"
        }
      ],
      "tag": "getmaximumiterationstotransport",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::getMaximumIterationsToTransport"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstotransport",
          "tag": "ogetmaximumiterationstotransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "iline0"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "iline0",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "iline0",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline0iextraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "iline0extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::iLine0Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstotransport",
          "tag": "ogetmaximumiterationstotransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "iline1"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "iline1",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "iline1",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline1iextraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "iline1extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::iLine1Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstotransport",
          "tag": "ogetmaximumiterationstotransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "iline2"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "iline2",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "iline2",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline2iextraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "iline2extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::iLine2Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstotransport",
          "tag": "ogetmaximumiterationstotransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "iline3"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "iline3",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "iline3",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline3iextraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "iline3extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::iLine3Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstotransport",
          "tag": "ogetmaximumiterationstotransport",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "png",
              "type": "FILE"
            }
          ],
          "tag": "oivisualization",
          "type": "OUTPUT"
        }
      ],
      "tag": "ivisualization",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::iVisualization"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "getmaximumiterationstotransport",
          "tag": "ogetmaximumiterationstotransport",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "t_step",
              "type": "NUMERIC"
            },
            {
              "name": "dt",
              "type": "NUMERIC"
            },
            {
              "name": "time",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "type": "NUMERIC"
            },
            {
              "name": "flow_l",
              "type": "NUMERIC"
            },
            {
              "name": "flow_n_linear_iterations",
              "type": "NUMERIC"
            },
            {
              "name": "flow_final_linear_residual",
              "type": "NUMERIC"
            },
            {
              "name": "flow_norm_delta",
              "type": "NUMERIC"
            },
            {
              "name": "flow_norm_delta_u",
              "type": "NUMERIC"
            },
            {
              "name": "flow_converged",
              "type": "TEXT"
            }
          ],
          "tag": "osolversimulationflow",
          "type": "OUTPUT"
        }
      ],
      "tag": "solversimulationflow",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::SolverSimulationFlow"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationflow",
          "tag": "osolversimulationflow",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "t_step",
              "type": "NUMERIC"
            },
            {
              "name": "dt",
              "type": "NUMERIC"
            },
            {
              "name": "time",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "type": "NUMERIC"
            },
            {
              "name": "transport_l",
              "type": "NUMERIC"
            },
            {
              "name": "transport_n_linear_iterations",
              "type": "NUMERIC"
            },
            {
              "name": "transport_final_linear_residual",
              "type": "NUMERIC"
            },
            {
              "name": "transport_norm_delta",
              "type": "NUMERIC"
            },
            {
              "name": "transport_norm_delta_u",
              "type": "NUMERIC"
            },
            {
              "name": "transport_converged",
              "type": "TEXT"
            }
          ],
          "tag": "osolversimulationtransport",
          "type": "OUTPUT"
        }
      ],
      "tag": "solversimulationtransport",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::SolverSimulationTransport"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "t_step",
              "type": "NUMERIC"
            },
            {
              "name": "time",
              "type": "NUMERIC"
            },
            {
              "name": "dt",
              "type": "NUMERIC"
            },
            {
              "name": "n_flow_linear_iterations_total",
              "type": "NUMERIC"
            },
            {
              "name": "n_flow_nonlinear_iterations_total",
              "type": "NUMERIC"
            },
            {
              "name": "n_transport_linear_iterations_total",
              "type": "NUMERIC"
            },
            {
              "name": "n_transport_nonlinear_iterations_total",
              "type": "NUMERIC"
            },
            {
              "name": "solution_converged",
              "type": "TEXT"
            },
            {
              "name": "error",
              "type": "NUMERIC"
            }
          ],
          "tag": "ocomputesolutionchange",
          "type": "OUTPUT"
        }
      ],
      "tag": "computesolutionchange",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::EvaluateTimeStepControl"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "computesolutionchange",
          "tag": "ocomputesolutionchange",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "t_step",
              "type": "NUMERIC"
            },
            {
              "name": "time",
              "type": "NUMERIC"
            },
            {
              "name": "dt",
              "type": "NUMERIC"
            },
            {
              "name": "ts_converged",
              "type": "TEXT"
            }
          ],
          "tag": "ocomputetimestep",
          "type": "OUTPUT"
        }
      ],
      "tag": "computetimestep",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::EvaluateTimeStepControl"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "first_step_refinement",
              "type": "TEXT"
            },
            {
              "name": "t_step",
              "type": "NUMERIC"
            },
            {
              "name": "before_n_active_elem",
              "type": "NUMERIC"
            },
            {
              "name": "after_n_active_elem",
              "type": "NUMERIC"
            }
          ],
          "tag": "omeshrefinement",
          "type": "OUTPUT"
        }
      ],
      "tag": "meshrefinement",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::MeshRefinement"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            }
          ],
          "tag": "omeshwriter",
          "type": "OUTPUT"
        }
      ],
      "tag": "meshwriter",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::MeshWriter"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "line0"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "line0",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "line0",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline0extraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "line0extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::Line0Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "line1"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "line1",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "line1",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline1extraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "line1extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::Line1Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "line2"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "line2",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "line2",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline2extraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "line2extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::Line2Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "extractors": [
            {
              "cartridge": "EXTRACTION",
              "extension": "CSV",
              "tag": "line3"
            }
          ],
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "u",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "v",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "w",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "p",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "s",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "d",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "r",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "vtkvalidpointmask",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "arc_length",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "points0",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "points1",
              "extractor": "line3",
              "type": "NUMERIC"
            },
            {
              "name": "points2",
              "extractor": "line3",
              "type": "NUMERIC"
            }
          ],
          "tag": "oline3extraction",
          "type": "OUTPUT"
        }
      ],
      "tag": "line3extraction",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::Line3Extraction"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "solversimulationtransport",
          "tag": "osolversimulationtransport",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "time_step",
              "type": "NUMERIC"
            },
            {
              "name": "png",
              "type": "FILE"
            }
          ],
          "tag": "ovisualization",
          "type": "OUTPUT"
        }
      ],
      "tag": "visualization",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::Visualization"
        }
      ]
    },
    {
      "sets": [
        {
          "dependency": "meshwriter",
          "tag": "omeshwriter",
          "type": "INPUT"
        },
        {
          "attributes": [
            {
              "name": "simulationid",
              "type": "NUMERIC"
            },
            {
              "name": "xdmf",
              "type": "FILE"
            },
            {
              "name": "n_processors",
              "type": "NUMERIC"
            }
          ],
          "tag": "omeshaggregator",
          "type": "OUTPUT"
        }
      ],
      "tag": "meshaggregator",
      "programs": [
        {
          "path": "/scratch/10061a/vitorss/simulation/sedimentation-rest/libmesh-sedimentation",
          "name": "libmesh-sedimentation-opt::MeshAggregator"
        }
      ]
    }
  ],
  "tag": "sedimentation"
}'