// Implement member functions for MeanCurvatureFlow class.
#include "Mem-3dg.h"

#include "Energy_Handler.h"
#include "Beads.h"
#include "Interaction.h"

#include "BeadGeometry.h"

#include <fstream>
#include <omp.h>
#include <sys/stat.h>

#include <chrono>
#include <Eigen/Core>
#include "geometrycentral/surface/remeshing.h"
// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>
/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
using namespace std;
Mem3DG::Mem3DG(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo)
{

  // Build member variables: mesh, geometry
  mesh = inputMesh;
  geometry = inputGeo;
  // velocity = VertexData<Vector3>(*mesh,{0,0,0});
  Old_norm2 = 0;
  Current_norm2 = 0;
  H_Vector_0 = VertexData<double>(*mesh, 0.0);
  dH_Vector = VertexData<double>(*mesh, 0.0);
  // H_target = VertexData<double> (*mesh,0.0);
  system_time = 0;
  grad_norm = 0.0;
  is_test = false;
  pulling = false;
  Area_evol_steps = 100000;
  stop_increasing = false;
  small_TS = false;
  recentering = true;
  boundary = false;
  Field = "None";
  Field_vals.resize(0);

  momentum = false;
  learn_rate = 0.5;

  BFGS_iter = 0;
  N_data = 0;
  mean_E = 0;
  var_E = 0;
  mean_Grad = 0;
  var_Grad = 0;
  Turn_normal_iter = false;
}
Mem3DG::Mem3DG(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, Bead input_Bead)
{

  // Build member variables: mesh, geometry
  mesh = inputMesh;
  geometry = inputGeo;
  // velocity = VertexData<Vector3>(*mesh,{0,0,0});
  Old_norm2 = 0;
  Current_norm2 = 0;
  H_Vector_0 = VertexData<double>(*mesh, 0.0);
  dH_Vector = VertexData<double>(*mesh, 0.0);
  // H_target = VertexData<double> (*mesh,0.0);
  system_time = 0;
  grad_norm = 0.0;
  is_test = false;
  Bead_1 = input_Bead;
  // Beads.push_back(input_Bead);
  pulling = false;
  Save_SS = false;
  stop_increasing = false;
  small_TS = false;
  remesh_flag = false;
}

Mem3DG::Mem3DG(ManifoldSurfaceMesh *inputMesh, VertexPositionGeometry *inputGeo, bool test)
{

  // Build member variables: mesh, geometry
  mesh = inputMesh;
  geometry = inputGeo;
  // velocity = VertexData<Vector3>(*mesh,{0,0,0});
  Old_norm2 = 0;
  Current_norm2 = 0;
  H_Vector_0 = VertexData<double>(*mesh, 0.0);
  dH_Vector = VertexData<double>(*mesh, 0.0);
  // H_target = VertexData<double> (*mesh,0.0);
  system_time = 0;
  grad_norm = 0.0;
  is_test = test;
  pulling = false;
  stop_increasing = false;
  small_TS = false;
  remesh_flag = false;
}

// Destructor: clear internal containers. This class does not own
// `mesh`, `geometry`, or `Sim_handler` so those pointers are not
// deleted here.
Mem3DG::~Mem3DG()
{
  Beads.clear();
  s_list.clear();
  y_list.clear();
  rho_list.clear();
}

void Mem3DG::Add_bead(Bead *bead)
{

  Beads.push_back(bead);
  Beads[Beads.size() - 1]->Reasign_mesh(mesh, geometry);
  // std::cout<<"Rpoi"
  // Beads[Beads.size()-1]->Bead_I = bead->Bead_I;
  // std::cout<<"Moving pointers\n";
  // std::cout<<"The directino of the I is "<< Beads[Beads.size()-1]->Bead_I << "MEDG \n";
}


// Recentering between steps. The membrane is moved back to the origin; with a
// field (and field_aware) its leftmost vertex is moved to x = -1 instead.
// Returns the value the beads are shifted back by: the old centre of mass, or
// the field displacement.
Vector3 Mem3DG::recenter_membrane(bool field_aware)
{
  if (field_aware && Field != "None")
  {
    Vector3 leftmost = Vector3({1e10, 1e10, 1e10});
    for (Vertex v : mesh->vertices())
    {
      if (geometry->inputVertexPositions[v].x < leftmost.x)
        leftmost = geometry->inputVertexPositions[v];
    }
    leftmost = Vector3({-1.0, 0.0, 0.0}) - leftmost;
    VertexData<Vector3> Displacement(*mesh, leftmost);
    geometry->inputVertexPositions += Displacement;
    return leftmost;
  }
  Vector3 CoM = geometry->centerOfMass();
  geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
  return CoM;
}

VertexData<Vector3> Mem3DG::SurfaceGrad() const
{

  size_t index;
  Vector3 Normal;
  size_t N_vert = mesh->nVertices();
  VertexData<Vector3> Force(*mesh);
  Vector3 u;
  Halfedge he_grad;

  for (Vertex v : mesh->vertices())
  {
    Normal = {0, 0, 0};
    Force[v] = {0, 0, 0};
    Force[v] = -1 * 2 * geometry->vertexNormalMeanCurvature(v);
  }

  // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
  return Force;
}


/**
 * @brief Performs the backtracking algorithm for the Mem3DG class.
 *
 *
 * @return The optimal step size for the backtracking algorithm.
 */
double Mem3DG::Backtracking()
{
  // std::cout << "Backtracking\n";

  double c1 = 0.5;
  double rho = 0.5;
  double alpha = 1;
  double position_Projeection = 0;
  double X_pos;

  // We calculate the dot product between the 2 gradients
  if (momentum && system_time > 0 && Sim_handler->Previous_grad.size() > 0)
  {
    double product = 0;
    for (Vertex v : mesh->vertices())
    {
      product += dot(Sim_handler->Current_grad[v], Sim_handler->Previous_grad[v]);
    }
    if (product > 0)
    {
      Sim_handler->Current_grad = Sim_handler->Current_grad + Sim_handler->Previous_grad * learn_rate;
      for (Bead *b : Beads)
        b->Total_force = b->Total_force + b->Prev_Total_force * learn_rate; // This line needs to be tested still
    }
  }
  double previousE = 0;
  Sim_handler->Calculate_energies(&previousE);
  for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
  {
    if (isnan(Sim_handler->Energy_values[i]))
      std::cout << "Energy " << Sim_handler->Energies[i] << " is nan\n";
  }

  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  if (recentering)
    recenter_membrane(true);
  Vector3 CoM = geometry->centerOfMass();

  initial_pos = geometry->inputVertexPositions;

  std::vector<Vector3> Bead_init;

  for (size_t i = 0; i < Beads.size(); i++)
    Bead_init.push_back(Beads[i]->Pos);

  double Projection = 0;
  Vector3 center;

  for (Vertex v : mesh->vertices())
  {
    if (isnan(Sim_handler->Current_grad[v].x) || isnan(Sim_handler->Current_grad[v].y) || isnan(Sim_handler->Current_grad[v].z))
      // std::cout << " Is this force nan at vertex but norm2 is " << Sim_handler->Current_grad[v.getIndex()].norm2() << "\n";
      std::cout << v.getIndex() << " ";
  }
  // std::cout << " \n";

  geometry->inputVertexPositions += alpha * Sim_handler->Current_grad;
  bool nanflag = false;
  for (Vertex v : mesh->vertices())
  {
    if (isnan(geometry->inputVertexPositions[v].norm2()))
      nanflag = true;
  }

  if (nanflag)
    std::cout << "At least one vertex has nan position\n";

  // geometry->refreshQuantities();
  center = geometry->centerOfMass();
  Vector3 Vertex_pos;

  // We move the beads;
  for (size_t i = 0; i < Beads.size(); i++)
    Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));

  Projection = Sim_handler->Gradient_norms[Sim_handler->Energies.size()];

  size_t bead_count = 0;
  NewE = 0.0;
  Sim_handler->Calculate_energies(&NewE);
  size_t counter = 0;
  bool displacement_cond = true;

  if (Projection < 1e-7)
  {
    small_TS = true;
    std::cout << "The energy diff is quite small and so is the gradient\n";
    std::cout << "The energy diff is" << abs(NewE - previousE) / previousE << "\n";
    std::cout << "The projection is" << Projection << "\n";
    return -1.0;
  }

  while (true)
  {
    displacement_cond = true;

    for (size_t i = 0; i < Beads.size(); i++)
      displacement_cond = displacement_cond && Beads[i]->Total_force.norm() * alpha < 0.1 * Beads[i]->sigma;

    if (NewE <= previousE - c1 * alpha * Projection && displacement_cond && fabs(NewE - previousE) < 1e2)
    {
      break;
    }

    if (std::isnan(NewE))
    {
      std::cout << "When backtracking newE is nan";
      alpha = -1.0;
      break;
    }

    alpha *= rho;

    if ((fabs((NewE - previousE) / previousE) < 1e-4 && Projection < 0.05) || Projection < 1e-5)
    {
      small_TS = true;
      std::cout << "The energy diff is quite small and so is the gradient\n";
      std::cout << "The energy diff is" << abs(NewE - previousE) << "\n";
      std::cout << "The projection is" << Projection << "\n";
      return -1.0;
    }

    if (alpha < 1e-10)
    {
      BFGS_iter = 0;
      remesh_flag = true;
      // std::cout << "The remesh flag is true\n";
      // std::cout << "THe timestep got small so the simulation would end \n";
      // std::cout << "THe timestep is " << alpha << " \n";
      // std::cout << "The energy diff is" << abs(NewE - previousE) << "\n";
      // std::cout << "THe relative energy diff  is" << abs((NewE - previousE) / previousE) << "\n";
      // std::cout << "The projection is" << Projection << "\n";
      // std::cout << "The projection is too big " << (Projection > 1.0e8) << " \n";
      if (Projection > 1.0e8)
      {
        // return alpha;
        std::cout << "The gradient got crazy\n";
        std::cout << "The projection is " << Projection << "\n";
        // geometry->inputVertexPositions = initial_pos;
        return -1;
      }
      small_TS = true;
      if (Projection >= 5000)
        small_TS = false;

      break;
    }
    else if (small_TS)
      small_TS = false;
    if (alpha > 0)
    {
      geometry->inputVertexPositions = initial_pos + alpha * Sim_handler->Current_grad;

      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
        Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));
      }
    }
    else
    {
      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
      }
      geometry->inputVertexPositions = initial_pos;
    }
    bead_count = 0;
    NewE = 0.0;
    Sim_handler->Calculate_energies(&NewE);
  }

  nanflag = false;

  for (Vertex v : mesh->vertices())
    if (isnan(geometry->inputVertexPositions[v].x + geometry->inputVertexPositions[v].y + geometry->inputVertexPositions[v].z))
      nanflag = true;

  if (nanflag)
    std::cout << "After backtracking one vertex is nan :( also the value of alpha is" << alpha << " \n";
  if (alpha <= 0.0)
  {
    geometry->inputVertexPositions = initial_pos;
  }
  if (recentering)
  {
    if (alpha <= 0.0)
      std::cout << "NotRecentering after crisis\n";
    else
    {
      CoM = recenter_membrane(true);
      for (size_t i = 0; i < Beads.size(); i++)
        Beads[i]->Pos -= CoM;
    }
  }

  return alpha;
}


// Line search of the Newton methods on the norm of the Lagrangian gradient.
// normal = true for the normal-direction variant (Newton-Normal).
double Mem3DG::Backtracking_newton(Eigen::VectorXd p_lambda, double Projection, double Current_grad_norm, bool normal)
{

  double c1 = 0.2;
  double rho = 0.5;
  double alpha = 1;
  // alpha = 5e-4;
  double position_Projeection = 0;
  double X_pos;

  int N_vert = mesh->nVertices();
  int N_beads = Beads.size();

  VertexData<Vector3> Step_newton = Sim_handler->Current_grad;
  std::vector<Vector3> Step_beads(0);

  for (int bi = 0; bi < N_beads; bi++)
    Step_beads.push_back(Beads[bi]->Total_force);

  double PrevNorm = Current_grad_norm;

  double NewNorm;

  VertexData<Vector3> initial_pos(*mesh);
  Eigen::VectorXd initial_lag;
  std::vector<Vector3> initial_bead_pos(0);

  // We will open a file to log the backtracking process
  std::ofstream backtrack_log;

  backtrack_log.open(basic_name + "backtrack_log.txt", std::ios::app);
  // So the order will be : Prevnorm Projection NEWNORM NEWNORM NEWNORM ...

  backtrack_log << PrevNorm << " ";
  if (!normal)
    backtrack_log << Projection << " ";

  if (recentering)
    recenter_membrane(true);
  Vector3 CoM = geometry->centerOfMass();

  // std::cout<<"Saving initial condition \n";

  initial_pos = geometry->inputVertexPositions;
  initial_lag = Sim_handler->Lagrange_mult;

  for (int bi = 0; bi < N_beads; bi++)
  {
    initial_bead_pos.push_back(Beads[bi]->Pos);
  }

  std::vector<Vector3> Bead_init;

  for (size_t i = 0; i < N_beads; i++)
    Bead_init.push_back(Beads[i]->Pos);

  // double Projection = 0;
  Vector3 center;

  // std::cout<<"stepping\n";

  geometry->inputVertexPositions += alpha * Step_newton;
  // We move the beads;
  for (size_t i = 0; i < Beads.size(); i++)
    Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));

  // for(int bi = 0; bi < N_beads; bi++) Beads

  // We only change the lagrange multipliers of the constraint of volume and area

  for (int i = 0; i < Sim_handler->N_constraints; i++)
  {
    if (Sim_handler->Constraints[i] == "Volume" || Sim_handler->Constraints[i] == "Area")
    {
      Sim_handler->Lagrange_mult[i] -= alpha * p_lambda[i];
    }
  }

  // std::cout<<"Done stepping\n";

  bool nanflag = false;
  for (Vertex v : mesh->vertices())
  {
    if (isnan(geometry->inputVertexPositions[v].norm2()))
      nanflag = true;
  }

  if (nanflag)
    std::cout << "At least one vertex has nan position\n";

  // geometry->refreshQuantities();
  center = geometry->centerOfMass();
  Vector3 Vertex_pos;

  size_t bead_count = 0;
  NewNorm = 0.0;
  NewNorm = lagrangian_norm(normal);

  size_t counter = 0;

  bool displacement_cond = true;

  // std::cout<<"THe projection is " << Projection <<" \n";
  // std::cout<<"The number of beads is " << Beads.size() << " \n";

  if (Projection < 1e-7)
  {
    small_TS = true;
    std::cout << "The norm  diff is quite small and so is the gradient\n";
    std::cout << "The norm diff is" << abs(NewNorm - PrevNorm) / PrevNorm << "\n";
    std::cout << "The projection is" << Projection << "\n";
    return -1.0;
  }

  while (true)
  {
    displacement_cond = true;
    // std::cout<<"The Norms are Prev: " << PrevNorm << " New: " << NewNorm << "\n";

    backtrack_log << NewNorm << " ";

    for (size_t i = 0; i < Beads.size(); i++)
      displacement_cond = displacement_cond && Beads[i]->Total_force.norm() * alpha < 0.1 * Beads[i]->sigma;

    // if( fabs(PrevNorm-NewNorm) <= alpha * Projection && NewNorm < PrevNorm  ) {
    if (NewNorm < PrevNorm)
    {

      if (fabs(NewNorm - PrevNorm) > 5e1 && false)
      {

        std::cout << "The energies are ";
        for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
          std::cout << Sim_handler->Energies[i] << " is " << Sim_handler->Energy_values[i] << " ";
        std::cout << " \n";
        std::cout << "The projection is " << Projection << " \n";

        double Max_projection = 0.0;
        int maxproj_index = 0;

        double maxDisplacement = 0.0;
        std::cout << "Finding breaking point\n";
        for (Vertex v : mesh->vertices())
        {
          if (Sim_handler->Current_grad[v].norm() > Max_projection)
          {
            Max_projection = Sim_handler->Current_grad[v].norm();
            maxproj_index = v.getIndex();
          }
          double displacement = (geometry->inputVertexPositions[v] - initial_pos[v]).norm();
          if (displacement > maxDisplacement)
          {
            maxDisplacement = displacement;
          }
        }
        std::cout << "The max displacement is " << maxDisplacement << " \n";
        std::cout << "The value of alpha is " << alpha << " \n";
        std::cout << "We will recalculate the energies, lets go back one step for now\n";
        geometry->inputVertexPositions = initial_pos;

        Sim_handler->Lagrange_mult = initial_lag;

        // geometry->refreshQuantities();
        mesh->compress();
        // I want something else

        alpha = 0.0;

        // Lets troubleshoot this hehe
        std::cout << "The previous energy was" << PrevNorm << " \n";
        std::cout << "The projection of the bigges vertex is " << Max_projection << " \n";
        std::cout << "This vertex is located at " << geometry->inputVertexPositions[maxproj_index] << " \n";
        std::cout << "This vertex in init pos is  at " << initial_pos[maxproj_index] << " \n";

        // Lets explore the sorroundings
        Vertex v = mesh->vertex(maxproj_index);
        for (Face f : v.adjacentFaces())
        {
          std::cout << "The adjacent faces are " << f.getIndex() << " \n";
          std::cout << "With area " << geometry->faceArea(f) << " \n";
        }
        for (Halfedge he : v.outgoingHalfedges())
        {
          std::cout << "The adjacent halfedges are " << he.getIndex() << " \n";
          std::cout << "With cotan " << geometry->cotan(he) << " \n";
        }
      }
      backtrack_log << "\n";
      backtrack_log.close();
      break;
    }

    if (std::isnan(NewNorm))
    {
      std::cout << "Grad norm is nan \n";
      backtrack_log << "\n";
      backtrack_log.close();
      alpha = -1.0;
      break;
    }

    alpha *= rho;
    if ((abs((NewNorm - PrevNorm) / PrevNorm) < 1e-7 && Projection < 0.5) || Projection < 1e-5)
    {
      small_TS = true;
      std::cout << "The energy diff is quite small and so is the gradient\n";
      std::cout << "The energy diff is" << abs(NewNorm - PrevNorm) / PrevNorm << "\n";
      std::cout << "The projection is" << Projection << "\n";
      backtrack_log << "\n";
      backtrack_log.close();
      return -1.0;
    }

    if (alpha < 1e-10)
    {
      // std::cout << "THe timestep got small so the simulation would end \n";
      // std::cout << "THe timestep is " << alpha << " \n";
      // std::cout << "The NORM diff is" << abs(NewNorm - PrevNorm) << "\n";
      // std::cout << "THe relative energy diff  is" << abs((NewNorm - PrevNorm) / PrevNorm) << "\n";
      // std::cout << "The projection is" << Projection << "\n";
      // std::cout << "The projection is too big " << (Projection > 1.0e8) << " \n";
      if (Projection > 1.0e8)
      {
        // return alpha;
        std::cout << "The gradient got crazy\n";
        std::cout << "The projections is " << Projection << "\n";
        geometry->inputVertexPositions = initial_pos;
        Sim_handler->Lagrange_mult = initial_lag;
        backtrack_log << "\n";
        backtrack_log.close();
        return -1;
      }
      if (Projection < 100)
      {
        small_TS = true;
      }

      break;

      // LEts try to step when it get super small
    }

    else if (small_TS)
      small_TS = false;
    // std::cout<<"System time is" << system_time <<" \n";
    if (alpha > 0)
    {
      // std::cout<<"UPDATING POSITIONS\n";
      geometry->inputVertexPositions = initial_pos + alpha * Step_newton;

      for (int i = 0; i < Sim_handler->N_constraints; i++)
      {
        if (Sim_handler->Constraints[i] == "Volume" || Sim_handler->Constraints[i] == "Area")
        {
          // std::cout<<"Updating lagrange multipliers\n";
          Sim_handler->Lagrange_mult[i] = initial_lag[i] - alpha * p_lambda[i];
        }
      }
      // std::cout<<"The lagrange multipliers are " << Sim_handler->Lagrange_mult.transpose() << "\n";

      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
        Beads[i]->Total_force = Step_beads[i];
        Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));
      }
    }
    else
    {
      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
        // Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
      }
      geometry->inputVertexPositions = initial_pos;
      Sim_handler->Lagrange_mult = initial_lag;
    }

    // // geometry->refreshQuantities();

    bead_count = 0;

    // NewE = 0.0;
    // Sim_handler->Calculate_energies(&NewE);
    NewNorm = 0.0;
    NewNorm = lagrangian_norm(normal);
  }

  backtrack_log << "\n";
  backtrack_log.close();

  nanflag = false;

  for (Vertex v : mesh->vertices())
    if (isnan(geometry->inputVertexPositions[v].x + geometry->inputVertexPositions[v].y + geometry->inputVertexPositions[v].z))
      nanflag = true;

  if (nanflag)
    std::cout << "After backtracking one vertex is nan :( also the value of alpha is" << alpha << " \n";
  if (alpha <= 0.0)
  {
    // std::cout<<"Repositioning\n";
    geometry->inputVertexPositions = initial_pos;
  }
  if (recentering)
  {
    if (alpha <= 0.0)
      std::cout << "NotRecentering after crisis\n";
    else
    {
      CoM = recenter_membrane(true);
      for (size_t i = 0; i < Beads.size(); i++)
        Beads[i]->Pos -= CoM;
    }
  }

  // std::cout<<"The difference in energy is " << fabs(NewE-previousE) <<"(: \n";
  // std::cout<<"The new norm is " << NewNorm << "\n";
  return alpha;
}

// The Calculate_Lag_norm functions accumulate into their argument
double Mem3DG::lagrangian_norm(bool normal)
{
  double norm = 0.0;
  if (normal)
  {
    Sim_handler->Calculate_Lag_norm_Normal(&norm);
    return norm;
  }
  Sim_handler->Calculate_Lag_norm(&norm);
  return 0.5 * norm;
}

double Mem3DG::Backtracking_grad(Eigen::VectorXd p_lambda, double Projection, double Current_grad_norm)
{
  return Backtracking_newton(p_lambda, Projection, Current_grad_norm, false);
}

double Mem3DG::Backtracking_grad_Normal(Eigen::VectorXd p_lambda, double Projection, double Current_grad_norm)
{
  return Backtracking_newton(p_lambda, Projection, Current_grad_norm, true);
}


/**
 * @brief Performs the backtracking algorithm for the Mem3DG class using a Force given.
 *
 * This function is responsible for performing the backtracking algorithm in the Mem3DG class. It takes in the force vector, a vector of energy names, and a matrix of energy constants as input. It updates the position of the vertices and beads based on the force, and calculates the energy values for each energy term. It then performs backtracking to find the optimal step size that minimizes the energy.
 *
 * @param Force The force vector applied to the vertices.
 * @param Energies A vector of energy names.
 * @param Energy_constants A matrix of energy constants.
 * @return The optimal step size for the backtracking algorithm.
 */
double Mem3DG::Backtracking_BFGS(VertexData<Vector3> Force, std::vector<Vector3> Bead_forces)
{

  double c1 = 1e-4;
  double rho = 0.5;
  double alpha = 1;
  double position_Projeection = 0;
  double X_pos;

  // std::cout << "Backtracking\n";
  double previousE = 0;
  Sim_handler->Calculate_energies(&previousE);

  for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
  {
    if (isnan(Sim_handler->Energy_values[i]))
      std::cout << "Energy " << Sim_handler->Energies[i] << " is nan\n";
  }

  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  if (recentering)
    recenter_membrane(false);
  Vector3 CoM = geometry->centerOfMass();

  initial_pos = geometry->inputVertexPositions;

  std::vector<Vector3> Bead_init;

  for (size_t i = 0; i < Beads.size(); i++)
    Bead_init.push_back(Beads[i]->Pos);

  double Projection = 0;
  Vector3 center;

  for (Vertex v : mesh->vertices())
  {
    if (isnan(Force[v].x) || isnan(Force[v].y) || isnan(Force[v].z))
    {
      std::cout << BFGS_iter << " " << v.getIndex() << " ";
    }
    // std::cout << " Is this force nan at vertex but norm2 is " << Force[v.getIndex()].norm2() << "\n";
  }
  // std::cout << "\n";

  // std::cout<<"doing the stepping \n";
  geometry->inputVertexPositions += alpha * Force;
  bool nanflag = false;
  for (Vertex v : mesh->vertices())
  {
    if (isnan(geometry->inputVertexPositions[v].norm2()))
    {
      // std::cout << "Nan vertex index " << v.getIndex() << " \n";
      // std::cout << "Boundary? " << v.isBoundary() << "\n";
      nanflag = true;
      break;
    }
  }

  if (nanflag)
  {
    nanflag = true;
    std::cout << "At least one vertex has nan position\n";
  }
  center = geometry->centerOfMass();
  Vector3 Vertex_pos;

  // We move the beads;
  // std::cout<<"Moving beads\n";
  for (size_t i = 0; i < Beads.size(); i++)
    Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}), Bead_forces[i]);
  Projection = 0.0;
  // std::cout<<"Calculating projection\n";
  for (Vertex v : mesh->vertices())
  {
    Projection += Force[v].norm2();
  }

  size_t bead_count = 0;
  NewE = 0.0;
  Sim_handler->Calculate_energies(&NewE);
  size_t counter = 0;

  bool displacement_cond = true;

  // std::cout<<"THe projection is " << Projection <<" \n";
  // std::cout<<"The number of beads is " << Beads.size() << " \n";

  if (Projection < 1e-5)
  {
    small_TS = true;
    std::cout << "The energy diff is quite small and so is the gradieeent\n";
    std::cout << "The energy diff is" << abs(NewE - previousE) / previousE << "\n";
    std::cout << "The projection is" << Projection << "\n";
    return -1.0;
  }
  // std::cout<<"Backtracking\n";
  while (true)
  {
    displacement_cond = true;

    for (size_t i = 0; i < Beads.size(); i++)
      displacement_cond = displacement_cond && Beads[i]->Total_force.norm() * alpha < 0.1 * Beads[i]->sigma;

    if (NewE <= previousE - c1 * alpha * Projection && displacement_cond && fabs(NewE - previousE) < 1e2)
    {

      break;
    }

    if (std::isnan(NewE))
    {
      std::cout << "The new Energy is Nan\n";
      alpha = -1.0;
      break;
    }

    alpha *= rho;
    if ((abs((NewE - previousE) / previousE) < 1e-6 && Projection < 1e-6) || Projection < 1e-6)
    {
      small_TS = true;
      std::cout << "The energy diff is quite small and so is the gradient\n";
      std::cout << "The energy diff is" << abs((NewE - previousE) / previousE) << "\n";
      std::cout << "The projection is" << Projection << "\n";
      std::cout << "Finishing sim\n";
      return -1.0;
    }

    if (alpha < 1e-10)
    {
      // BFGS_iter = 0;

      if (Projection > 1.0e5)
      {
        remesh_flag = true;
        // return alpha;
        std::cout << "The gradient got crazy\n";
        std::cout << "The projections is " << Projection << "\n";
        geometry->inputVertexPositions = initial_pos;
        BFGS_iter = -1;
        // return -1; //I may need to change this later
      }
      // small_TS = true;

      break;
    }
    else if (small_TS)
      small_TS = false;
    // std::cout<<"System time is" << system_time <<" \n";
    if (alpha > 0)
    {
      geometry->inputVertexPositions = initial_pos + alpha * Force;

      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
        Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}), Bead_forces[i]);
      }
    }
    else
    {
      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
        // Beads[i]->Move_bead(alpha, Vector3({0,0,0}));
      }
      geometry->inputVertexPositions = initial_pos;
    }

    // geometry->refreshQuantities();

    bead_count = 0;

    NewE = 0.0;
    Sim_handler->Calculate_energies(&NewE);
  }

  // std::cout<<"The energy final is" << NewE <<" \n";

  nanflag = false;

  for (Vertex v : mesh->vertices())
    if (isnan(geometry->inputVertexPositions[v].x + geometry->inputVertexPositions[v].y + geometry->inputVertexPositions[v].z))
      nanflag = true;

  if (nanflag)
  {
    std::cout << "After backtracking one vertex is nan :( also the value of alpha is" << alpha << " \n";
    std::cout << "The energy values are ";
    for (size_t i = 0; i < Sim_handler->Energy_values.size(); i++)
    {
      std::cout << Sim_handler->Energy_values[i] << " ";
    }
    std::cout << "\n";
  }

  if (alpha <= 0.0)
  {
    // std::cout<<"Repositioning\n";
    // geometry->inputVertexPositions = initial_pos;
  }
  if (recentering)
  {
    if (alpha <= 0.0)
      std::cout << "NotRecentering after crisis\n";
    else
    {
      CoM = recenter_membrane(false);
      for (size_t i = 0; i < Beads.size(); i++)
        Beads[i]->Pos -= CoM;
    }
  }

  //

  // std::cout<<"The difference in energy is " << fabs(NewE-previousE) <<"(: \n";

  return alpha;
}


// One row of Output_data.txt:
// time step Volume Area <energy values> Total_E <gradient norms> backtrackstep
void Mem3DG::write_output_row(std::ofstream &Sim_data, double t, double Volume, double Area, double tot_E, double step) const
{
  Sim_data << t << " " << discreteTs << " " << Volume << " " << Area << " ";
  for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
    Sim_data << Sim_handler->Energy_values[i] << " ";
  Sim_data << tot_E << " ";
  for (size_t i = 0; i < Sim_handler->Gradient_norms.size(); i++)
    Sim_data << Sim_handler->Gradient_norms[i] << " ";
  Sim_data << step << " \n";
}

// One row per bead in Bead_<i>_data.txt: step position force
void Mem3DG::write_bead_rows(const std::vector<std::string> &Bead_data_filenames) const
{
  for (size_t i = 0; i < Beads.size(); i++)
  {
    std::ofstream Bead_data(Bead_data_filenames[i], std::ios_base::app);
    Bead_data << discreteTs << " " << Beads[i]->Pos.x << " " << Beads[i]->Pos.y << " " << Beads[i]->Pos.z << " " << Beads[i]->Total_force.x << " " << Beads[i]->Total_force.y << " " << Beads[i]->Total_force.z << " \n";
  }
}

double Mem3DG::integrate(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data)
{

  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  // std::cout<<"Got to integrating\n";
  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;

  A = 0;
  V = 0;

  start = chrono::steady_clock::now();
  Sim_handler->Calculate_gradient();
  end = chrono::steady_clock::now();

  time_gradients = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  double alpha = 1e-3;
  double backtrackstep;
  double Grad_tot_norm = 0;
  double r_eff = 0;

  start = chrono::steady_clock::now();
  // std::cout<<"Backtracking?\n";
  if (backtrack)
  {
    backtrackstep = Backtracking();
  }
  else
  {
    double some_E;
    Sim_handler->Calculate_energies(&some_E);
    backtrackstep = timestep;
    geometry->inputVertexPositions += backtrackstep * Sim_handler->Current_grad;

    for (size_t i = 0; i < Beads.size(); i++)
    {
      Beads[i]->Move_bead(backtrackstep, Vector3({0, 0, 0}));
    }
  }

  end = chrono::steady_clock::now();
  time_backtracking = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  geometry->refreshQuantities();
  // Here we will refres the quantities. So the area calculation is a lil bit less expensive
  if (Save_output_data || backtrackstep < 0)
  {
    geometry->requireFaceAreas();
    A = geometry->totalArea();
    geometry->unrequireFaceAreas();
    V = geometry->totalVolume();
    double tot_E = 0;
    for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
      tot_E += Sim_handler->Energy_values[i];
    write_output_row(Sim_data, time, V, A, tot_E, backtrackstep);
  }

  if (Save_output_data)
    write_bead_rows(Bead_data_filenames);


  return backtrackstep;
}

// L-BFGS two-loop recursion: applies the inverse Hessian estimate stored in
// s_list/y_list/rho_list to q. q is overwritten by the first loop.
Eigen::VectorXd Mem3DG::lbfgs_direction(Eigen::VectorXd &q) const
{
  std::vector<double> alpha_list(m);
  for (int i = BFGS_iter - 1; i >= 0 && BFGS_iter - i < m; i--)
  {
    double alpha_i = rho_list[i % m] * s_list[i % m].dot(q);
    alpha_list[i % m] = alpha_i;
    q = q - alpha_i * y_list[i % m];
  }
  Eigen::VectorXd r = q;
  for (int i = std::max(0, BFGS_iter - m); i < BFGS_iter; i++)
  {
    double beta_i = rho_list[i % m] * y_list[i % m].dot(r);
    r = r + s_list[i % m] * (alpha_list[i % m] - beta_i);
  }
  return r;
}

double Mem3DG::integrate_BFGS_Normal(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data)
{
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  double backtrackstep;

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;

  size_t bead_count = 0;

  start = chrono::steady_clock::now();
  std::vector<Vector3> Bead_forces(Beads.size());
  if (BFGS_iter == 0)
  {
    Sim_handler->Calculate_gradient();

    VertexData<double> Grad_E = VertexData<double>(*mesh, 0.0);
    for (Vertex v : mesh->vertices())
    {
      Grad_E[v] = dot(Sim_handler->Current_grad[v], Sim_handler->Vertex_normals[v]);
    }
    VertexData<Vector3> Force(*mesh, Vector3{0.0, 0.0, 0.0});
    for (Vertex v : mesh->vertices())
    {
      Force[v] = Grad_E[v] * Sim_handler->Vertex_normals[v];
    }
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Bead_forces[bi] = Vector3{0.0, 0.0, 0.0};
    }

    backtrackstep = Backtracking_BFGS(Force, Bead_forces);

    Eigen::VectorXd Aux_vector = Eigen::VectorXd::Zero((mesh->nVertices()));

    for (Vertex v : mesh->vertices())
    {
      Aux_vector[v.getIndex()] = backtrackstep * Grad_E[v];
    }

    s_list.resize(0);
    y_list.resize(0);
    rho_list.resize(0);

    s_list.push_back(Aux_vector);
    // here the vertices are alreay updated no?

    geometry->refreshQuantities();

    Sim_handler->Calculate_gradient();

    for (Vertex v : mesh->vertices())
    {
      Aux_vector[v.getIndex()] = dot(Sim_handler->Current_grad[v], Sim_handler->Vertex_normals[v]) - dot(Sim_handler->Previous_grad[v], Sim_handler->Vertex_normals[v]);
    }
    y_list.push_back(Aux_vector);

    double rho_i = 1.0 / (s_list[0].dot(y_list[0]));
    if (isnan(rho_i) || isinf(rho_i) || fabs(s_list[0].dot(y_list[0])) < 1e-10)
    {
      std::cout << "\tThe value of rho_ i  is something " << rho_i << "\n";
      BFGS_iter = -1;
    }
    rho_list.push_back(rho_i);
  }
  else
  {
    // FOr some reason i am not calculating the gradient here?
    Eigen::VectorXd Grad_E = Eigen::VectorXd::Zero(mesh->nVertices());

    for (Vertex v : mesh->vertices())
    {
      Grad_E[v.getIndex()] = dot(Sim_handler->Current_grad[v], Sim_handler->Vertex_normals[v]);
    }

    Eigen::VectorXd r = lbfgs_direction(Grad_E);
    VertexData<Vector3> Force(*mesh, Vector3({0.0, 0.0, 0.0}));
    for (Vertex v : mesh->vertices())
    {
      Force[v] = r[v.getIndex()] * Sim_handler->Vertex_normals[v];
    }
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Bead_forces[bi] = Vector3{0.0, 0.0, 0.0};
    }

    backtrackstep = Backtracking_BFGS(Force, Bead_forces);

    Eigen::VectorXd s_k = backtrackstep * r;

    // geometry->refreshQuantities();
    Sim_handler->Calculate_gradient();
    for (Vertex v : mesh->vertices())
    {
      Grad_E[v.getIndex()] = dot((Sim_handler->Current_grad[v] - Sim_handler->Previous_grad[v]), Sim_handler->Vertex_normals[v]);
    }

    if (BFGS_iter < m)
    {
      s_list.push_back(s_k);
      y_list.push_back(Grad_E);
      rho_list.push_back(1.0 / (s_k.dot(Grad_E)));
    }
    else
    {
      s_list[BFGS_iter % m] = s_k;
      y_list[BFGS_iter % m] = Grad_E;

      rho_list[BFGS_iter % m] = 1.0 / (s_k.dot(Grad_E));
    }

    if (fabs(s_k.dot(Grad_E)) < 1e-10 || isinf(rho_list[BFGS_iter % m]) || isnan(rho_list[BFGS_iter % m]))
    {
      BFGS_iter = -1;
    }
  }
  // std::cout<<"DOne with iteration\n";
  BFGS_iter += 1;

  // HERE I CAN FREE sk AND gRAD_E

  double r_eff = 0;

  end = chrono::steady_clock::now();
  time_backtracking = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  double tot_E = 0;
  Sim_handler->Calculate_energies(&tot_E);

  N_data += 1;
  mean_E = mean_E + (tot_E - mean_E) / N_data;
  if (N_data == 1)
  {
    var_E = 0;
  }
  else
  {
    var_E = var_E * (N_data - 1) / (N_data) + (tot_E - mean_E) * (tot_E - mean_E) / (N_data - 1);
  }

  if (N_data % 200 == 0)
  {

    // In this case we will print the value of this quantity
    std::cout << "After" << N_data << "steps of iterations the mean Energy is " << mean_E << " and the variance is " << var_E << "\n";
    N_data = 0;
    mean_E = 0;
    var_E = 0;
  }

  if (Save_output_data || backtrackstep < 0.0)
  {
    V = geometry->totalVolume();
    A = geometry->totalArea();
    write_output_row(Sim_data, time, V, A, tot_E, backtrackstep);
  }

  if (Bead_data_filenames.size() != 0 && Save_output_data)
    write_bead_rows(Bead_data_filenames);

  return backtrackstep;
}

double Mem3DG::integrate_BFGS(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data)
{

  // We need to store the Hessian somewhere, maybe it should be received as a pointer
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  double backtrackstep;

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;
  Turn_normal_iter = false;

  size_t bead_count = 0;

  start = chrono::steady_clock::now();

  std::vector<Vector3> Bead_forces(Beads.size());

  if (BFGS_iter == 0)
  {
    Sim_handler->Calculate_gradient();
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Bead_forces[bi] = Beads[bi]->Total_force;
    }
    backtrackstep = Backtracking_BFGS(Sim_handler->Current_grad, Bead_forces);

    Eigen::VectorXd Aux_vector = Eigen::VectorXd::Zero(3 * (mesh->nVertices() + Beads.size()));
    for (Vertex v : mesh->vertices())
    {
      Aux_vector[3 * v.getIndex()] = backtrackstep * Sim_handler->Current_grad[v].x;
      Aux_vector[3 * v.getIndex() + 1] = backtrackstep * Sim_handler->Current_grad[v].y;
      Aux_vector[3 * v.getIndex() + 2] = backtrackstep * Sim_handler->Current_grad[v].z;
    }
    int N_vert = mesh->nVertices();
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Aux_vector[3 * N_vert + 3 * bi] = backtrackstep * Beads[bi]->Total_force.x;
      Aux_vector[3 * N_vert + 3 * bi + 1] = backtrackstep * Beads[bi]->Total_force.y;
      Aux_vector[3 * N_vert + 3 * bi + 2] = backtrackstep * Beads[bi]->Total_force.z;
    }
    s_list.resize(0);
    y_list.resize(0);
    rho_list.resize(0);

    s_list.push_back(Aux_vector);

    geometry->refreshQuantities();
    Sim_handler->Calculate_gradient();

    Vector3 Diff;
    for (Vertex v : mesh->vertices())
    {
      Diff = Sim_handler->Current_grad[v] - Sim_handler->Previous_grad[v];
      Aux_vector[3 * v.getIndex()] = Diff.x;
      Aux_vector[3 * v.getIndex() + 1] = Diff.y;
      Aux_vector[3 * v.getIndex() + 2] = Diff.z;
    }
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Diff = Beads[bi]->Total_force - Beads[bi]->Prev_Total_force;
      Aux_vector[3 * N_vert + 3 * bi] = Diff.x;
      Aux_vector[3 * N_vert + 3 * bi + 1] = Diff.y;
      Aux_vector[3 * N_vert + 3 * bi + 2] = Diff.z;
    }
    y_list.push_back(Aux_vector);

    double rho_i = 1.0 / (s_list[0].dot(y_list[0]));
    if (isnan(rho_i) || isinf(rho_i || fabs(s_list[0].dot(y_list[0])) < 1e-10))
    {
      std::cout << "\tThe value of rho_ i  is something " << rho_i << "\n";
      BFGS_iter = -1;
    }
    rho_list.push_back(rho_i);
  }
  else
  {

    Eigen::VectorXd Grad_vec = Eigen::VectorXd::Zero(3 * (mesh->nVertices() + Beads.size()));
    for (Vertex v : mesh->vertices())
    {
      Grad_vec[3 * v.getIndex()] = Sim_handler->Current_grad[v].x;
      Grad_vec[3 * v.getIndex() + 1] = Sim_handler->Current_grad[v].y;
      Grad_vec[3 * v.getIndex() + 2] = Sim_handler->Current_grad[v].z;
    }
    int N_vert = mesh->nVertices();
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Grad_vec[3 * N_vert + 3 * bi] = Beads[bi]->Total_force.x;
      Grad_vec[3 * N_vert + 3 * bi + 1] = Beads[bi]->Total_force.y;
      Grad_vec[3 * N_vert + 3 * bi + 2] = Beads[bi]->Total_force.z;
    }

    Eigen::VectorXd r = lbfgs_direction(Grad_vec);

    // So here we have r which is the product of the inverse Hessian with the gradient, we just need to do the backtracking and then update the lists
    VertexData<Vector3> Force(*mesh, Vector3({0.0, 0.0, 0.0}));
    for (Vertex v : mesh->vertices())
    {
      Force[v].x = r[3 * v.getIndex()];
      Force[v].y = r[3 * v.getIndex() + 1];
      Force[v].z = r[3 * v.getIndex() + 2];
    }
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
      Bead_forces[bi] = Vector3{r[3 * N_vert + 3 * bi], r[3 * N_vert + 3 * bi + 1], r[3 * N_vert + 3 * bi + 2]};
    }

    backtrackstep = Backtracking_BFGS(Force, Bead_forces);

    // Ok here we have a lil problem when.
    if (backtrackstep > 0.0)
    {
      Eigen::VectorXd s_k = backtrackstep * r;

      geometry->refreshQuantities();
      mesh->compress();
      Sim_handler->Calculate_gradient();
      for (Vertex v : mesh->vertices())
      {
        Grad_vec[3 * v.getIndex()] = Sim_handler->Current_grad[v].x - Sim_handler->Previous_grad[v].x;
        Grad_vec[3 * v.getIndex() + 1] = Sim_handler->Current_grad[v].y - Sim_handler->Previous_grad[v].y;
        Grad_vec[3 * v.getIndex() + 2] = Sim_handler->Current_grad[v].z - Sim_handler->Previous_grad[v].z;
      }
      for (size_t bi = 0; bi < Beads.size(); bi++)
      {
        Grad_vec[3 * N_vert + 3 * bi] = Beads[bi]->Total_force.x - Beads[bi]->Prev_Total_force.x;
        Grad_vec[3 * N_vert + 3 * bi + 1] = Beads[bi]->Total_force.y - Beads[bi]->Prev_Total_force.y;
        Grad_vec[3 * N_vert + 3 * bi + 2] = Beads[bi]->Total_force.z - Beads[bi]->Prev_Total_force.z;
      }

      if (BFGS_iter < m)
      {
        s_list.push_back(s_k);
        y_list.push_back(Grad_vec);
        rho_list.push_back(1.0 / (s_k.dot(Grad_vec)));
      }
      else
      {
        // In case this is not the case we do
        s_list[BFGS_iter % m] = s_k;
        y_list[BFGS_iter % m] = Grad_vec;
        rho_list[BFGS_iter % m] = 1.0 / (s_k.dot(Grad_vec));
      }
      if (fabs(s_k.dot(Grad_vec)) < 1e-10 || isinf(rho_list[BFGS_iter % m]))
      {
        // std::cout << "THe value of s_k dot is " << s_k.dot(Grad_vec) << " \n";
        // std::cout << "Resetting bfgs iters\n";
        BFGS_iter = -1;
      }
    }
  }

  BFGS_iter += 1;
  double r_eff = 0;

  end = chrono::steady_clock::now();
  time_backtracking = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  double tot_E = 0;
  Sim_handler->Calculate_energies(&tot_E);

  N_data += 1;
  mean_E = mean_E + (tot_E - mean_E) / N_data;

  // var_E = var_E*(N_data-1)/(N_data)  + (tot_E-mean_E)*(tot_E-mean_E)/(N_data-1);

  // double grad_norm = Sim_handler->Gradient_norms[Sim_handler->Gradient_norms.size()-1];
  // mean_Grad = mean_Grad + (grad_norm - mean_Grad)/N_data;
  if (N_data == 1)
  {
    var_E = 0;
  }
  else
  {
    var_E = var_E * (N_data - 1) / (N_data) + (tot_E - mean_E) * (tot_E - mean_E) / (N_data - 1);
    // var_Grad = var_Grad*(N_data-1)/(N_data)  + (grad_norm - mean_Grad)*(grad_norm - mean_Grad)/(N_data-1);
  }

  if (N_data % 100 == 0)
  {

    // // In this case we will print the value of this quantity
    std::cout << "After" << N_data << "steps of iterations the mean Energy is " << mean_E << " and the variance is " << var_E << "\n";
    // std::cout<<"After" << N_data << "steps of iterations the mean Energy is " << mean_Grad <<" and the variance is " << var_Grad << "\n";

    if (var_E < 0.0006)
      Turn_normal_iter = true;
    N_data = 0;
    mean_E = 0;
    var_E = 0;
  }

  if (Save_output_data || backtrackstep <= 0.0)
  {
    V = geometry->totalVolume();
    A = geometry->totalArea();
    write_output_row(Sim_data, time, V, A, tot_E, backtrackstep);
  }
  if (Bead_data_filenames.size() != 0 && (Save_output_data || backtrackstep < 0.0))
    write_bead_rows(Bead_data_filenames);
  return backtrackstep;
}

// Right-hand side of a constraint row in the Newton systems: minus the
// constraint violation (the centre of mass and rotation rows are zero).
double Mem3DG::constraint_rhs(const std::string &name, double area) const
{
  if (name == "Volume")
    return -1 * (geometry->totalVolume() - Sim_handler->Trgt_vol);
  if (name == "Area")
    return -1 * (area - Sim_handler->Trgt_area);
  return 0.0;
}

double Mem3DG::integrate_Newton(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames)
{

  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  bool regularize = true;
  // std::cout<<"Got to integrating\n";

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;

  size_t bead_count = 0;
  A = 0;
  geometry->requireFaceAreas();
  for (Face f : mesh->faces())
    A += geometry->faceAreas[f];

  // std::cout<<"AREA CALCULATED\n";

  size_t N_vert = mesh->nVertices();
  int N_beads = Beads.size();
  int N_constraints = 0;

  for (size_t i = 0; i < Constraints.size(); i++)
  {
    if (Constraints[i] == "Volume")
      N_constraints += 1;
    if (Constraints[i] == "Area")
      N_constraints += 1;
    if (Constraints[i] == "CMx")
      N_constraints += 1;
    if (Constraints[i] == "CMy")
      N_constraints += 1;
    if (Constraints[i] == "CMz")
      N_constraints += 1;
    if (Constraints[i] == "Rx")
      N_constraints += 1;
    if (Constraints[i] == "Ry")
      N_constraints += 1;
    if (Constraints[i] == "Rz")
      N_constraints += 1;
  }

  Sim_handler->N_constraints = N_constraints;
  // SparseMatrix<double> H2_mat;//(N_vert*3+N_constraints,N_vert*3+N_constraints);
  // SparseMatrix<double> H1_mat;
  SparseMatrix<double> Hessian;
  SparseMatrix<double> Hessian_constraints;

  Eigen::SimplicialLDLT<SparseMatrix<double>> solverHess;

  start = chrono::steady_clock::now();
  Sim_handler->Calculate_gradient();
  Sim_handler->Calculate_Jacobian();
  Hessian = Sim_handler->Calculate_Hessian();

  SparseMatrix<double> LHS(3 * (N_vert + N_beads) + N_constraints, 3 * (N_vert + N_beads) + N_constraints);
  Eigen::VectorXd RHS(3 * (N_vert + N_beads) + N_constraints);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int row = 0; row < Sim_handler->Jacobian_constraints.rows(); row++)
  {
    for (int col = 0; col < Sim_handler->Jacobian_constraints.cols(); col++)
    {
      // Ok so now we need to add this
      tripletList.push_back(T(col, row + 3 * (N_vert + N_beads), Sim_handler->Jacobian_constraints(row, col)));
      tripletList.push_back(T(row + 3 * (N_vert + N_beads), col, Sim_handler->Jacobian_constraints(row, col)));
    }
  }

  int row;
  int col;
  double value;
  int maxrow = 0;
  int maxcol = 0;

  int smol_counter = 0;
  for (long int k = 0; k < Hessian.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(Hessian, k); it; ++it)
    {
      value = it.value();
      row = it.row();
      col = it.col();
      if (value < 1e-8 && value > -1e-8)
        smol_counter += 1;
      tripletList.push_back(T(row, col, value));
    }
  }
  // We regularize the Hessian by adding a diagonal
  if (regularize)
  {
    for (size_t i = 0; i < 3 * (N_vert + N_beads); i++)
    {
      tripletList.push_back(T(i, i, 1e-6));
    }
  }

  Eigen::VectorXd LambdaJ = Sim_handler->Jacobian_constraints.transpose() * Sim_handler->Lagrange_mult;

  Vector3 Force;
  double dual_area = 0.0;
  for (size_t vi = 0; vi < mesh->nVertices(); vi++)
  {
    Force = Sim_handler->Current_grad[vi];
    RHS(3 * vi) = LambdaJ(3 * vi) + Force.x;
    RHS(3 * vi + 1) = LambdaJ(3 * vi + 1) + Force.y;
    RHS(3 * vi + 2) = LambdaJ(3 * vi + 2) + Force.z;
  }
  // I need to add the force of the bead here
  for (size_t bi = 0; bi < Sim_handler->Beads.size(); bi++)
  {
    Force = Sim_handler->Beads[bi]->Total_force;
    RHS(3 * (N_vert + bi)) = LambdaJ(3 * (N_vert + bi)) + Force.x;
    RHS(3 * (N_vert + bi) + 1) = LambdaJ(3 * (N_vert + bi) + 1) + Force.y;
    RHS(3 * (N_vert + bi) + 2) = LambdaJ(3 * (N_vert + bi) + 2) + Force.z;
  }

  for (int Ci = 0; Ci < N_constraints; Ci++)
    RHS(3 * (N_vert + N_beads) + Ci) = constraint_rhs(Constraints[Ci], A);

  LHS.setFromTriplets(tripletList.begin(), tripletList.end());
  solverHess.compute(LHS);

  Eigen::VectorXd result = solverHess.solve(RHS);

  // std::cout<<"Solved\n";
  VertexData<Vector3> Force_result(*mesh, Vector3({0.0, 0.0, 0.0}));
  // Eigen::VectorXd Grad_L = LambdaJ+ ;
  bool flag = false;
  for (size_t vi = 0; vi < mesh->nVertices(); vi++)
  {
    if (mesh->vertex(vi).isBoundary())
    {
      Force_result[vi] = Vector3({0.0, 0.0, 0.0});
      result(3 * vi) = 0.0;
      result(3 * vi + 1) = 0.0;
      result(3 * vi + 2) = 0.0;
      continue;
    }
    for (Vertex vj : mesh->vertex(vi).adjacentVertices())
    {
      if (vj.isBoundary())
      {
        Force_result[vi] = Vector3({0.0, 0.0, 0.0});
        result(3 * vi) = 0.0;
        result(3 * vi + 1) = 0.0;
        result(3 * vi + 2) = 0.0;
        flag = true;
        break;
      }
    }
    if (flag)
    {
      flag = false;
      continue;
    }
    Force_result[vi].x = result(3 * vi);
    Force_result[vi].y = result(3 * vi + 1);
    Force_result[vi].z = result(3 * vi + 2);
  }

  Sim_handler->Previous_grad = Sim_handler->Current_grad;
  Sim_handler->Current_grad = Force_result;

  for (int bi = 0; bi < N_beads; bi++)
  {
    Force = Vector3{result(3 * (N_vert + bi)), result(3 * (N_vert + bi) + 1), result(3 * (N_vert + bi) + 2)};
    Beads[bi]->Total_force = Force;
  }

  double Projection = result.transpose() * LHS * RHS;
  double Current_grad_norm = 0.5 * RHS.dot(RHS);

  double backtrackstep;
  if (false)
  {
    // if(result.dot(RHS)<0 || Projection <0.0){
    //
    std::cout << "The result is not a descent direction, we will not backtrack\n";
    small_TS = true;
    backtrackstep = 0.0;
    // backtrackstep = integrate(Sim_data,time,Bead_data_filenames,Save_output_data);
  }
  else
  {

    if (backtrack)
    {
      backtrackstep = Backtracking_grad(result.tail(N_constraints), Projection, Current_grad_norm);
    }
    else
      backtrackstep = timestep;

    // if(backtrackstep < 2e-6) small_TS = true;
    geometry->refreshQuantities();
    double TotE = 0.0;

    Sim_handler->Calculate_energies(&TotE);
    A = 0;
    geometry->requireFaceAreas();
    for (Face f : mesh->faces())
      A += geometry->faceAreas[f];
    geometry->unrequireFaceAreas();

    write_output_row(Sim_data, time + backtrackstep, geometry->totalVolume(), A, TotE, backtrackstep);
  }

  if (Bead_data_filenames.size() != 0 && Save_output_data)
    write_bead_rows(Bead_data_filenames);

  return backtrackstep;
}

VertexData<Vector3> Mem3DG::Newton_Normal_step(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames)
{
  std::cout << "THe lagrange multipliers are " << Sim_handler->Lagrange_mult.transpose() << "\n";
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  bool regularize = true;

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;

  size_t bead_count = 0;
  geometry->requireFaceAreas();
  A = 0.0;
  for (Face f : mesh->faces())
    A += geometry->faceAreas[f];
  geometry->unrequireFaceAreas();

  size_t N_vert = mesh->nVertices();
  int N_beads = Beads.size();
  int N_constraints = 0;

  for (size_t i = 0; i < Constraints.size(); i++)
  {
    if (Constraints[i] == "Volume")
      N_constraints += 1;
    if (Constraints[i] == "Area")
      N_constraints += 1;
    if (Constraints[i] == "CMx")
      N_constraints += 1;
    if (Constraints[i] == "CMy")
      N_constraints += 1;
    if (Constraints[i] == "CMz")
      N_constraints += 1;
    if (Constraints[i] == "Rx")
      N_constraints += 1;
    if (Constraints[i] == "Ry")
      N_constraints += 1;
    if (Constraints[i] == "Rz")
      N_constraints += 1;
  }
  std::cout << "The number of constraints is " << N_constraints << "\n";
  Sim_handler->N_constraints = N_constraints;
  Sim_handler->Constraints = Constraints;
  // SparseMatrix<double> H2_mat;//(N_vert*3+N_constraints,N_vert*3+N_constraints);
  // SparseMatrix<double> H1_mat;
  SparseMatrix<double> Hessian;
  SparseMatrix<double> Hessian_constraints;

  Eigen::SimplicialLDLT<SparseMatrix<double>> solverHess;

  start = chrono::steady_clock::now();
  // std::cout << "Calculating gradients\n";
  Sim_handler->Calculate_gradient();
  Hessian = Sim_handler->Calculate_Hessian_Normal();
  Sim_handler->Calculate_Jacobian_Normal();

  SparseMatrix<double> LHS((N_vert) + N_constraints, (N_vert) + N_constraints);
  Eigen::VectorXd RHS((N_vert) + N_constraints);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int row = 0; row < Sim_handler->Jacobian_constraints.rows(); row++)
  {
    for (int col = 0; col < Sim_handler->Jacobian_constraints.cols(); col++)
    {
      tripletList.push_back(T(col, row + (N_vert), Sim_handler->Jacobian_constraints(row, col)));
      tripletList.push_back(T(row + (N_vert), col, Sim_handler->Jacobian_constraints(row, col)));
    }
  }

  int row;
  int col;
  double value;

  int maxrow = 0;
  int maxcol = 0;

  int smol_counter = 0;
  for (long int k = 0; k < Hessian.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(Hessian, k); it; ++it)
    {
      value = it.value();
      row = it.row();
      col = it.col();
      if (value < 1e-8 && value > -1e-8)
        smol_counter += 1;
      tripletList.push_back(T(row, col, value));
    }
  }
  // We regularize the Hessian by adding a diagonal
  if (regularize)
  {
    for (size_t i = 0; i < (N_vert + N_beads); i++)
    {
      tripletList.push_back(T(i, i, 1e-6));
    }
  }
  Eigen::VectorXd LambdaJ = Sim_handler->Jacobian_constraints.transpose() * Sim_handler->Lagrange_mult;

  std::cout << "The lagrange multipliers are " << Sim_handler->Lagrange_mult.transpose() << "\n";
  std::cout << "The size of Jacobian constraints is " << Sim_handler->Jacobian_constraints.rows() << " and " << Sim_handler->Jacobian_constraints.cols() << "\n";

  std::cout << "THe size of Lambda J is " << LambdaJ.rows() << " and " << LambdaJ.cols() << "\n";
  std::cout << "THe maximum value of Lambda J is " << LambdaJ.maxCoeff() << " and the minimum is " << LambdaJ.minCoeff() << "\n";
  Vector3 Force;
  double dual_area = 0.0;
  for (size_t vi = 0; vi < mesh->nVertices(); vi++)
  {
    Force = Sim_handler->Current_grad[vi];
    RHS(vi) = LambdaJ(vi) + dot(Force, Sim_handler->Vertex_normals[vi]);
  }
  std::cout << "RHS force from vertices ready\n";
  std::cout << "The size of RHS is " << RHS.rows() << " and " << RHS.cols() << "\n";

  for (int Ci = 0; Ci < N_constraints; Ci++)
    RHS((N_vert) + Ci) = constraint_rhs(Constraints[Ci], A);
  std::cout << "THe maximum value of RHS is " << RHS.maxCoeff() << " and the minimum is " << RHS.minCoeff() << "\n";

  LHS.setFromTriplets(tripletList.begin(), tripletList.end());
  solverHess.compute(LHS);

  Eigen::VectorXd result = solverHess.solve(RHS);

  std::cout << "The size of result is " << result.rows() << " and " << result.cols() << "\n";
  std::cout << "THe max value of result is " << result.maxCoeff() << " and the minimum is " << result.minCoeff() << "\n";

  VertexData<Vector3> Force_result(*mesh, Vector3({0.0, 0.0, 0.0}));
  // Eigen::VectorXd Grad_L = LambdaJ+ ;
  bool flag = false;
  for (size_t vi = 0; vi < mesh->nVertices(); vi++)
  {
    if (mesh->vertex(vi).isBoundary())
    {
      Force_result[vi] = Vector3({0.0, 0.0, 0.0});
      result(vi) = 0.0;
      continue;
    }
    for (Vertex vj : mesh->vertex(vi).adjacentVertices())
    {
      if (vj.isBoundary())
      {
        Force_result[vi] = Vector3({0.0, 0.0, 0.0});
        result(vi) = 0.0;
        flag = true;
        break;
      }
    }
    if (flag)
    {
      flag = false;
      continue;
    }
    Force_result[vi] = result(vi) * Sim_handler->Vertex_normals[vi];
  }

  for (int bi = 0; bi < N_beads; bi++)
  {
    Force = Vector3{0.0, 0.0, result((N_vert + bi))};

    Beads[bi]->Total_force = Force;
  }

  return Force_result;
  //
}
double Mem3DG::integrate_Newton_Normal(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames)
{
  std::cout << "THe lagrange multipliers are " << Sim_handler->Lagrange_mult.transpose() << "\n";
  auto start = chrono::steady_clock::now();
  auto end = chrono::steady_clock::now();
  auto construction_start = chrono::steady_clock::now();
  auto construction_end = chrono::steady_clock::now();
  auto solve_start = chrono::steady_clock::now();
  auto solve_end = chrono::steady_clock::now();

  bool regularize = true;

  double time_construct = 0;
  double time_solve = 0;
  double time_compute = 0;
  double time_gradients = 0;
  double time_backtracking = 0;

  size_t bead_count = 0;
  geometry->requireFaceAreas();
  A = 0.0;
  for (Face f : mesh->faces())
    A += geometry->faceAreas[f];
  geometry->unrequireFaceAreas();

  size_t N_vert = mesh->nVertices();
  int N_beads = Beads.size();
  int N_constraints = 0;

  for (size_t i = 0; i < Constraints.size(); i++)
  {
    if (Constraints[i] == "Volume")
      N_constraints += 1;
    if (Constraints[i] == "Area")
      N_constraints += 1;
    if (Constraints[i] == "CMx")
      N_constraints += 1;
    if (Constraints[i] == "CMy")
      N_constraints += 1;
    if (Constraints[i] == "CMz")
      N_constraints += 1;
    if (Constraints[i] == "Rx")
      N_constraints += 1;
    if (Constraints[i] == "Ry")
      N_constraints += 1;
    if (Constraints[i] == "Rz")
      N_constraints += 1;
  }
  std::cout << "The number of constraints is " << N_constraints << "\n";
  Sim_handler->N_constraints = N_constraints;
  Sim_handler->Constraints = Constraints;
  // SparseMatrix<double> H2_mat;//(N_vert*3+N_constraints,N_vert*3+N_constraints);
  // SparseMatrix<double> H1_mat;
  SparseMatrix<double> Hessian;
  SparseMatrix<double> Hessian_constraints;

  Eigen::SimplicialLDLT<SparseMatrix<double>> solverHess;

  start = chrono::steady_clock::now();
  // std::cout << "Calculating gradients\n";
  Sim_handler->Calculate_gradient();
  // Hessian = Sim_handler->Calculate_Hessian_Normal();
  Hessian = Sim_handler->Calculate_Hessian_Normal_clipped();
  Sim_handler->Calculate_Jacobian_Normal();

  SparseMatrix<double> LHS((N_vert) + N_constraints, (N_vert) + N_constraints);
  Eigen::VectorXd RHS((N_vert) + N_constraints);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int row = 0; row < Sim_handler->Jacobian_constraints.rows(); row++)
  {
    for (int col = 0; col < Sim_handler->Jacobian_constraints.cols(); col++)
    {
      tripletList.push_back(T(col, row + (N_vert), Sim_handler->Jacobian_constraints(row, col)));
      tripletList.push_back(T(row + (N_vert), col, Sim_handler->Jacobian_constraints(row, col)));
    }
  }

  int row;
  int col;
  double value;

  int maxrow = 0;
  int maxcol = 0;

  int smol_counter = 0;
  for (long int k = 0; k < Hessian.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(Hessian, k); it; ++it)
    {
      value = it.value();
      row = it.row();
      col = it.col();
      if (value < 1e-8 && value > -1e-8)
        smol_counter += 1;
      tripletList.push_back(T(row, col, value));
    }
  }
  // We regularize the Hessian by adding a diagonal
  if (regularize)
  {
    double safety_shift = 1e-10 * Hessian.diagonal().cwiseAbs().maxCoeff();
    std::cout << "Safety shift if " << safety_shift << " \n";
    for (size_t i = 0; i < (N_vert + N_beads); i++)
    {
      tripletList.push_back(T(i, i, safety_shift));
    }
  }
  Eigen::VectorXd LambdaJ = Sim_handler->Jacobian_constraints.transpose() * Sim_handler->Lagrange_mult;

  std::cout << "The lagrange multipliers are " << Sim_handler->Lagrange_mult.transpose() << "\n";
  std::cout << "The size of Jacobian constraints is " << Sim_handler->Jacobian_constraints.rows() << " and " << Sim_handler->Jacobian_constraints.cols() << "\n";

  std::cout << "THe size of Lambda J is " << LambdaJ.rows() << " and " << LambdaJ.cols() << "\n";
  std::cout << "THe maximum value of Lambda J is " << LambdaJ.maxCoeff() << " and the minimum is " << LambdaJ.minCoeff() << "\n";
  Vector3 Force;
  double dual_area = 0.0;
  for (size_t vi = 0; vi < mesh->nVertices(); vi++)
  {
    Force = Sim_handler->Current_grad[vi];
    RHS(vi) = LambdaJ(vi) + dot(Force, Sim_handler->Vertex_normals[vi]);
  }
  std::cout << "RHS force from vertices ready\n";
  std::cout << "The size of RHS is " << RHS.rows() << " and " << RHS.cols() << "\n";

  for (int Ci = 0; Ci < N_constraints; Ci++)
    RHS((N_vert) + Ci) = constraint_rhs(Constraints[Ci], A);
  std::cout << "THe maximum value of RHS is " << RHS.maxCoeff() << " and the minimum is " << RHS.minCoeff() << "\n";

  LHS.setFromTriplets(tripletList.begin(), tripletList.end());
  solverHess.compute(LHS);
  std::cout << "Is the solve a success?\n";
  if (solverHess.info() != Eigen::Success)
  {
    std::cout << "THe solve is not a succes?\n";
  }

  Eigen::VectorXd result = solverHess.solve(RHS);

  std::cout << "The size of result is " << result.rows() << " and " << result.cols() << "\n";
  std::cout << "THe max value of result is " << result.maxCoeff() << " and the minimum is " << result.minCoeff() << "\n";

  VertexData<Vector3> Force_result(*mesh, Vector3({0.0, 0.0, 0.0}));
  // Eigen::VectorXd Grad_L = LambdaJ+ ;
  bool flag = false;
  for (size_t vi = 0; vi < mesh->nVertices(); vi++)
  {
    if (mesh->vertex(vi).isBoundary())
    {
      Force_result[vi] = Vector3({0.0, 0.0, 0.0});
      result(vi) = 0.0;
      continue;
    }
    for (Vertex vj : mesh->vertex(vi).adjacentVertices())
    {
      if (vj.isBoundary())
      {
        Force_result[vi] = Vector3({0.0, 0.0, 0.0});
        result(vi) = 0.0;
        flag = true;
        break;
      }
    }
    if (flag)
    {
      flag = false;
      continue;
    }
    Force_result[vi] = result(vi) * Sim_handler->Vertex_normals[vi];
  }

  for (int bi = 0; bi < N_beads; bi++)
  {
    Force = Vector3{0.0, 0.0, result((N_vert + bi))};

    Beads[bi]->Total_force = Force;
  }

  //

  double Projection = result.transpose() * LHS * RHS;
  double Current_grad_norm = 0.0;

  Sim_handler->Calculate_Lag_norm_Normal(&Current_grad_norm);

  Sim_handler->Current_grad = Force_result;

  // Current grad norm is

  double backtrackstep;
  // if (false)
  // {

  // Ok so here
  double dirDerivative = result.head(mesh->nVertices()).dot(RHS.head(mesh->nVertices()));
  std::cout << "The directional derivative is " << dirDerivative << " \n";

  if (result.dot(RHS) < 0 && false)
  {
    //
    std::cout << "The result is not a descent direction, we will not backtrack\n";
    small_TS = true;
    backtrackstep = 0.0;
    // backtrackstep = integrate(Sim_data,time,Bead_data_filenames,Save_output_data);
  }
  else
  {

    if (backtrack)
    {
      // std::cout << "Projection is " << Projection << " and the grad norm is " << Current_grad_norm << "\n";
      backtrackstep = Backtracking_grad_Normal(result.tail(N_constraints), Projection, Current_grad_norm);
    }
    else
    {
      backtrackstep = timestep;
      geometry->inputVertexPositions += Force_result * backtrackstep;
    }
    std::cout << "THe backtrackstep is" << backtrackstep << " \n";
    geometry->refreshQuantities();
    double TotE = 0.0;

    Sim_handler->Calculate_energies(&TotE);
    geometry->requireFaceAreas();
    A = geometry->totalArea();
    write_output_row(Sim_data, time + backtrackstep, geometry->totalVolume(), A, TotE, backtrackstep);
  }

  if (Bead_data_filenames.size() != 0 && Save_output_data)
    write_bead_rows(Bead_data_filenames);

  return backtrackstep;
}


// The exponent will be 1e-6


VertexData<Vector3> Mem3DG::Grad_Bead(std::ofstream &Gradient_file, bool Save, bool Projection)
{

  std::cout << "The interaction to consider is " << Bead_1.interaction << "\n";
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  VertexData<Vector3> Finite_grad(*mesh);
  initial_pos = geometry->inputVertexPositions;

  // VertexData<Vector3> Gradients(*mesh,0.0);
  // double V=geometry->totalVolume();
  double E_bead = 0.0;
  double E_bead_back = 0.0;
  double E_bead_front = 0.0;
  double dr;
  double r_dist;
  // double D_P=-P0*(V-V_bar)/V_bar/V_bar;
  double total_grad_finite = 0;
  double total_grad_theory = 0;

  Vector3 grad{0.0, 0.0, 0.0};
  Vector3 grad_theory;
  Vector3 difference;
  Vector3 Area_grad;
  Vector3 r;
  E_bead = Bead_1.Energy();

  // E_vol=E_Pressure(D_P,V,V_bar);

  VertexData<Vector3> Calc_grad = Bead_1.Gradient();
  VertexData<Vector3> Grad_area = SurfaceGrad();
  dr = 1e-7;

  size_t N_vert = mesh->nVertices();
  std::cout << "The number of vertices is " << N_vert << " \n";
  // for(size_t index=0; index<N_vert; index++){
  size_t index;
  for (Vertex v : mesh->vertices())
  {
    // std::cout<<"One vertex\n";
    index = v.getIndex();
    grad_theory = Calc_grad[v];

    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    E_bead_front = Bead_1.Energy();

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    E_bead_back = Bead_1.Energy();

    grad.x = (E_bead_front - E_bead_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    E_bead_back = Bead_1.Energy();

    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    E_bead_front = Bead_1.Energy();

    grad.y = (E_bead_front - E_bead_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    E_bead_back = Bead_1.Energy();

    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    E_bead_front = Bead_1.Energy();

    grad.z = (E_bead_front - E_bead_back) / (2 * dr);

    if (Save)
    {
      difference = grad + grad_theory;


      Gradient_file << difference.x << " " << difference.y << " " << difference.z << " " << difference.norm() / grad.norm() << " " << difference.norm() << " " << grad.norm() / grad_theory.norm() << " \n"; //<< dot(r.unit(),difference.unit())<<" " << dot(Area_grad,difference.unit())<<" "<< dot(Area_grad,r.unit())<<" \n";//<< dot(r,HN) <<" \n" ;
      // difference= grad_theory;
      // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
      total_grad_theory += grad_theory.norm2();
      total_grad_finite += grad.norm2();
    }
    Finite_grad[v] = -1 * grad;
    geometry->inputVertexPositions[v] = initial_pos[v];
    // geometry->refreshQuantities();
  }
  if (Save)
  {
    Gradient_file << sqrt(total_grad_theory) << " " << sqrt(total_grad_finite) << " " << sqrt(total_grad_theory / total_grad_finite) << " \n";
  }

  return Finite_grad;
}

VertexData<Vector3> Mem3DG::Grad_tot_Area(std::ofstream &Gradient_file, bool Save) const
{
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  VertexData<Vector3> Finite_grad(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // VertexData<Vector3> Gradients(*mesh,0.0);
  double A = geometry->totalArea();

  double E_area = 0.0;
  double E_tot_area = 0.0;
  double E_area_back = 0.0;
  double E_area_front = 0.0;
  double dr;
  // double lambda=KA*(A-A_bar )/A_bar;
  double total_grad_finite = 0;
  double total_grad_theory = 0;

  Vector3 grad{0.0, 0.0, 0.0};
  Vector3 grad_theory;
  Vector3 difference;

  E_area = A;

  VertexData<Vector3> Calc_grad = SurfaceGrad();

  dr = 1e-7;

  size_t N_vert = mesh->nVertices();
  // for(size_t index=0; index<N_vert; index++){
  size_t index;
  for (Vertex v : mesh->vertices())
  {
    // A=geometry->totalArea();

    // E_tot_area=E_Surface(KA,A,A_bar);
    // std::cout<<"THe difference in energy is "<<E_tot_area-E_area<<" \n";
    index = v.getIndex();
    grad_theory = Calc_grad[v];
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    E_area_front = A;
    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    E_area_back = A;
    grad.x = (E_area_front - E_area_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();

    E_area_back = A;
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    E_area_front = A;
    grad.y = (E_area_front - E_area_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    E_area_back = A;
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    E_area_front = A;
    grad.z = (E_area_front - E_area_back) / (2 * dr);
    if (Save)
    {
      difference = grad + grad_theory;
      Gradient_file << difference.x / grad.norm() << " " << difference.y / grad.norm() << " " << difference.z / grad.norm() << " " << difference.norm() / grad.norm() << " " << grad.norm() / grad_theory.norm() << " \n";
      // difference= grad_theory;
      // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
      total_grad_theory += grad_theory.norm2();
      total_grad_finite += grad.norm2();
    }
    Finite_grad[v] = -1 * grad;
    geometry->inputVertexPositions[v] = initial_pos[v];
    // geometry->refreshQuantities();
  }
  if (Save)
  {
    Gradient_file << sqrt(total_grad_theory) << " " << sqrt(total_grad_finite) << " " << sqrt(total_grad_theory / total_grad_finite) << " \n";
  }

  return Finite_grad;
}


Eigen::Matrix2d Mem3DG::Face_sizing(Face f)
{
  Eigen::Matrix2d Sizing; //({0.0,0.0,0.0,0.0});
  // Eigen::Matrix2d Sizing{{0.0,0.0,0.0,0.0}};
  Sizing << 0.0, 0.0, 0.0, 0.0;

  Eigen::Matrix3d Rotation;

  Vector3 Normal1 = geometry->faceNormal(f);
  Vector3 Axis1({0.0, 1.0, 0.0});

  // std::cout<<"THe norm of the normal is " << Normal1.norm2()<<" \n";
  Eigen::Vector3d Normal{Normal1.x, Normal1.y, Normal1.z};
  Eigen::Vector3d Axis{Axis1.x, Axis1.y, Axis1.z};

  Eigen::Vector3d uvw = Normal.cross(Axis);

  double rcos = Normal.dot(Axis);
  double rsin = uvw.norm();

  if (rsin > 1e-7)
  {
    uvw /= rsin;
  }

  // Eigen::Matrix3d V_x{{0, uvw(2), -1*uvw(1) },{-1*uvw(2), 0, uvw(0)},{uvw(1), -1*uvw(0), 0}};

  // Eigen::Matrix3d V_x{{0, -1*uvw(2), uvw(1) },{uvw(2), 0, -1*uvw(0)},{-1*uvw(1), uvw(0), 0}};
  // Eigen::Matrix3d V_x({0, -1*uvw(2), uvw(1) ,uvw(2), 0, -1*uvw(0),-1*uvw(1), uvw(0), 0});
  Eigen::Matrix3d V_x;

  V_x << 0, -1 * uvw(2), uvw(1), uvw(2), 0, -1 * uvw(0), -1 * uvw(1), uvw(0), 0;

  Rotation = rcos * Eigen::Matrix3d::Identity() + rsin * V_x + uvw * uvw.transpose() * (1 - rcos);

  // We want to compare the outer proudct of two vectors
  // Eigen::Matrix3d outer{{uvw(0)*uvw(0),uvw(1)*uvw(0),uvw(2)*uvw(0) }, {uvw(0)*uvw(1),uvw(1)*uvw(1),uvw(2)*uvw(1)}, {uvw(0)*uvw(2),uvw(1)*uvw(2),uvw(2)*uvw(2)}};

  // std::cout<<"Outer is \n";
  // std::cout<<outer <<"\n";
  // std::cout<<"And the eigen one is "<<" \n";
  // std::cout<< uvw*uvw.transpose() <<"\n";

  Eigen::Vector3d trial({1.0, 5.0, 4.0});

  // std::cout<<"v_x is " << V_x << " \n";

  // std::cout<<"Trial has norm  " <<trial.norm() <<" \n";

  trial = Rotation * trial;
  //
  // std::cout<<"Trial has now norm " << trial.norm() << "\n ";
  // std::cout<<"Trial now is " << trial <<" \n";

  Vector3 Pos;
  Eigen::Vector3d spare_vec;
  vector<Eigen::Vector2d> Positions(0);
  vector<int> index(0);

  for (Vertex v : f.adjacentVertices())
  {

    Pos = geometry->inputVertexPositions[v];
    spare_vec = Eigen::Vector3d{Pos.x, Pos.y, Pos.z};
    spare_vec = Rotation * spare_vec;
    index.push_back(v.getIndex());
    Positions.push_back({spare_vec(0), spare_vec(2)});
    // std::cout<<"THe y pos is " << spare_vec(1) <<"\t";
  }
  // std::cout<<"\n";

  // std::cout<<"The order is " << index[0] << " then " << index[1] << "then " << index[2]<<" \n";

  int counter = 0;
  for (Edge e : f.adjacentEdges())
  {

    // std::cout<<"The original edgelength is " << geometry->edgeLength(e) << " \n";

    Eigen::Vector2d e_mat = Positions[(counter + 1) % 3] - Positions[counter % 3];

    // std::cout<<"THe projected length is " << e_mat.norm() << " \n";

    Eigen::Vector2d t_mat{-1 * e_mat(1), e_mat(0)};
    t_mat.normalize();

    double theta = geometry->dihedralAngle(e.halfedge());
    if (e.halfedge().face() != f)
    {
      theta = geometry->dihedralAngle(e.halfedge().twin());
    }

    Sizing -= 1 / 2. * theta * e_mat.norm() * t_mat * t_mat.transpose();

    counter += 1;
  }

  Sizing /= geometry->faceArea(f);

  return Sizing;
}

FaceData<double> Mem3DG::Face_sizings()
{

  FaceData<double> sizing(*mesh, 0.0);
  double aspect_min = 0.2;
  double refine_angle = 0.7;

  // Lets create the face sizing
  Eigen::Matrix2d Face_sizing_mat; //({0.0,0.0, 0.0,0.0});
  Face_sizing_mat << 0.0, 0.0, 0.0, 0.0;
  Eigen::EigenSolver<Eigen::Matrix2d> Solve;
  Eigen::Vector2d eigenvals;

  for (Face f : mesh->faces())
  {
    // For every face i need to do the sizing thingy
    Face_sizing_mat = Face_sizing(f);

    Face_sizing_mat = Face_sizing_mat.transpose() * Face_sizing_mat / (refine_angle * refine_angle);

    Solve.compute(Face_sizing_mat);

    eigenvals = Solve.eigenvalues().real();

    //
    for (int i = 0; i < 2; i++)
      eigenvals[i] = clamp(eigenvals[i], 1.f / (0.2 * 0.2), 1.f / (0.001 * 0.001));

    double lmax = std::max(eigenvals[0], eigenvals[1]);
    double lmin = lmax * aspect_min * aspect_min;

    for (int i = 0; i < 2; i++)
      if (eigenvals[i] < lmin)
        eigenvals[i] = lmin;

    sizing[f.getIndex()] = std::max(fabs(eigenvals[0]), fabs(eigenvals[1]));
  }

  return sizing;
}

VertexData<double> Mem3DG::Vert_sizing(FaceData<double> Face_sizings)
{

  VertexData<double> sizing(*mesh, 0.0);

  // We need to create the vertsizing first.
  // WE have the sizing of every face
  for (Vertex v : mesh->vertices())
  {
    double sum = 0.0;
    for (Face f : v.adjacentFaces())
    {
      sum += Face_sizings[f] * geometry->faceArea(f) / 3.;
    }
    sizing[v] = sum / geometry->vertexDualArea(v);
  }

  return sizing;
}


bool Mem3DG::Area_sanity_check()
{
  bool all_positive = true;
  for (Vertex v : mesh->vertices())
  {
    if (geometry->circumcentricDualArea(v) < 0)
    {
      all_positive = false;
    }
  }

  return all_positive;
}
