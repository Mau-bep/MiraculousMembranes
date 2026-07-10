// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>

// #include <omp.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>

#include <dirent.h>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/remeshing.h"

#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

#include <chrono>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

#include "args/args.hxx"
#include "imgui.h"

#include "Mem-3dg.h"
#include "Beads.h"
#include "Energy_Handler.h"
#include "math.h"

#include "io.hpp"
#include "simulation.hpp"
#include "conf.hpp"
#include "log.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace std;

template <class T, class... Args>
std::unique_ptr<T> make_unique(Args &&...args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh *mesh;
VertexPositionGeometry *geometry;

std::vector<std::unique_ptr<Interaction>> Interaction_container;

polyscope::PointCloud *psCloud;
SimplePolygonMesh Saved_mesh;

int counter_saved = 0;
std::vector<Vector3> Saved_beadpos(6);
VertexData<Vector3> Saved_vertex_positions;

// Polyscope visualization handle, to quickly add data to the surface

// polyscope::SurfaceMesh* psMesh;

// Some global variables
float TIMESTEP = -4;
float kappa = 1.0;

bool first_newton = true;

float H0 = 1.0;
float V_bar = (4 / 3) * PI * 10;
float A_bar;
double dA;
double Area;
double origA;
bool Area_constraint = false;
size_t minStepsReq = 0;

int current_t = 0;
int integration_counter = 0;
int integration_steps = 10;
int remesh_every = 1;
int last_remesh = 0;
// int last_save = 0;
float angleThresh = 0.2;
int save_interval;
bool integrating = false;
bool updateArea = true;
bool continueLastCall = false;
int savedSteps = 0;
int remesh_op_last = 0;
int trgt_remesh_op = 100;
double integral_error = 0;
bool adapt_remesh = true;
bool Count_remesh = false;
std::ofstream Remeshing_count;
int inspect_timestep = 0;
std::string Switch;
std::string filename;
size_t Switch_t = 0;

Eigen::VectorXd Lagrange_mults;
int remeshing_ops;

VertexData<Vector3> ORIG_VPOS; // original vertex positions
Vector3 CoM;                   // original center of mass
float Sx = 1;
float Sy = 1;
float Sz = 1;
float Tx = 0.0;
float Ty = 0.0;
float Tz = 0.0;

Mem3DG M3DG;
E_Handler Sim_handler;
VertexData<Vector3> NewtonStepSaved;
float NewtonTs = 1.0;

std::string basic_name;
int Save_slot = 0;
int Load_slot = 0;

std::vector<std::string> Energies(0);
std::vector<std::vector<double>> Energy_constants(0);
std::vector<std::vector<double>> Og_Energy_constants(0);
std::vector<double> Constants(0);
std::ofstream Sim_data;
std::vector<std::string> Bead_filenames;
std::ofstream Bead_datas;
// std::vector<Be>
std::vector<std::string> Constraints(0);
std::vector<Bead> Beads;
std::vector<std::string> bonds;
std::vector<std::vector<double>> constants;
// std::vector<Interaction *> Bead_Interactions(0);

std::vector<std::string> mesh_names(0);
std::vector<int> frames(0);
std::vector<std::vector<std::pair<size_t, Vector3>>> frameBeadPos;
std::vector<Vector3> CurrentBeadPos;
std::vector<Vector3> CurrentBeadForce;

polyscope::SurfaceMesh *psMesh;
double TOTAL_ANGLE_DEFECT;
size_t EULER_CHARACTERISTIC;

std::map<int, std::vector<std::array<double, 3>>> sColors;
std::map<int, std::vector<std::array<double, 3>>> sColors2;

polyscope::SurfaceVertexColorQuantity *vertexColors;
polyscope::SurfaceFaceColorQuantity *faceColors;
double vertexRadius;
double edgeRadius;

RemeshOptions Options;

enum shadeTpe
{
    sSHADED,
    sVertsizing,
    sFacesizing
};

std::vector<std::string> sNames = {"Off", "VertSizing", "FaceSizing"};

VertexData<double> Vert_sizings;
FaceData<double> Face_sizings;
EdgeData<double> Edge_sizings;

std::vector<Vector3> Gradient_vector;
VertexData<Vector3> Gradient_vertex;

std::vector<std::array<double, 3>> Vert_sizing_C;
std::vector<std::array<double, 3>> F_sizings_C;
std::vector<std::array<double, 3>> shaded_C;

std::string Integration = "Gradient_descent";
double scaling_factor = 1.0;

int frame_counter = 0;

std::vector<Eigen::VectorXd> Eigenvectors_Bending(0);
std::vector<Eigen::VectorXd> Eigenvectors_Surface(0);

std::array<double, 3> BLUE = {0.11, 0.388, 0.89};
glm::vec<3, float> ORANGE_VEC = {1, 0.65, 0};
std::array<double, 3> ORANGE = {1, 0.65, 0};

struct Deterministic_sort
{
    inline bool operator()(const std::pair<size_t, Vector3> &left, const std::pair<size_t, Vector3> &right)
    {
        return left.first < right.first;
    }
} deterministic_sort;

void showSelected()
{

    // Show selected vertices.
    geometry->requireVertexPositions();

    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it)
    {
        vertPos.push_back(geometry->vertexPositions[*it]);
    }
    polyscope::SurfaceGraphQuantity *showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor(ORANGE_VEC);

    // Show selected edges.
    // std::vector<Vector3> edgePos;
    // std::vector<std::array<size_t, 2>> edgeInd;
    // for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
    //      it != polyscope::state::subset.edges.end(); ++it)
    // {
    //     // Edge e = mesh->edge(*it);
    //     // std::cout << "Showing edge " << *it << "\n";
    //     // std::cout << "Connected by " << e.firstVertex() << " and " << e.secondVertex() << "\n";
    //     // std::cout<<"The vertices are" <<
    //     edgePos.push_back(geometry->vertexPositions[e.firstVertex()]);
    //     edgePos.push_back(geometry->vertexPositions[e.secondVertex()]);
    //     size_t i = edgeInd.size();
    //     edgeInd.push_back({2 * i, 2 * i + 1});
    // }
    // polyscope::SurfaceGraphQuantity *showEdges = psMesh->addSurfaceGraphQuantity("selected edges", edgePos, edgeInd);
    // showEdges->setEnabled(true);
    // showEdges->setRadius(edgeRadius);
    // showEdges->setColor(ORANGE_VEC);

    // Show selected faces.
    std::vector<std::array<double, 3>> faceColors(mesh->nFaces());
    for (size_t i = 0; i < mesh->nFaces(); i++)
    {
        faceColors[i] = BLUE;
    }
    for (std::set<size_t>::iterator it = polyscope::state::subset.faces.begin();
         it != polyscope::state::subset.faces.end(); ++it)
    {
        Face f = mesh->face(*it);
        std::cout << "The face " << *it << " has edges " << f.halfedge().edge().getIndex() << " " << f.halfedge().next().edge().getIndex() << " " << f.halfedge().next().next().edge().getIndex() << "\n";
        std::cout << "AND VERTICES " << f.halfedge().vertex().getIndex() << " " << f.halfedge().next().vertex().getIndex() << " " << f.halfedge().next().next().vertex().getIndex() << " \n";

        faceColors[*it] = ORANGE;
    }
    polyscope::SurfaceFaceColorQuantity *showFaces = psMesh->addFaceColorQuantity("selected faces", faceColors);
    showFaces->setEnabled(true);
    geometry->unrequireVertexPositions();
}

void redraw()
{
    geometry->requireVertexPositions();
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    vector<Vector3> BeadPositions(0);
    for (size_t i = 0; i < Beads.size(); i++)
    {
        BeadPositions.push_back(Beads[i].Pos);
        std::cout << "Positions is " << Beads[i].Pos << "\n";
    }
    if (Beads.size() > 0)
        psCloud->updatePointPositions(BeadPositions);
    polyscope::requestRedraw();
}

void Save_mesh(std::string basic_name, size_t current_t)
{
    Vector3 Pos;
    std::ofstream o(basic_name + "membrane_" + std::to_string(current_t) + ".obj");
    o << "#This is a meshfile from a saved state\n";
    for (Vertex v : mesh->vertices())
    {
        Pos = geometry->inputVertexPositions[v];
        o << "v " << Pos.x << " " << Pos.y << " " << Pos.z << "\n";
    }
    // Saving faces
    for (Face f : mesh->faces())
    {
        o << "f";

        for (Vertex v : f.adjacentVertices())
        {
            o << " " << v.getIndex() + 1;
        }
        o << "\n";
    }
    return;
}

void Save_dihedrals(std::string basic_name)
{
    // Build member variables: mesh, geometry
    Vector3 Pos;
    std::ofstream o(basic_name + "dihedrals_evol.txt", std::ios::app);
    for (Edge e : mesh->edges())
    {
        o << geometry->dihedralAngle(e.halfedge()) << " ";
    }
    o << "\n";
    o.close();
    return;
}

void Save_edgelengths(std::string basic_name)
{
    std::ofstream o(basic_name + "edgelengths.txt", std::ios::app);

    for (Edge e : mesh->edges())
    {
        o << geometry->edgeLength(e) << " ";
    }
    o << "\n";
    o.close();
    return;
}

SimplePolygonMesh save_geometry(ManifoldSurfaceMesh *mesh, VertexPositionGeometry *geometry)
{

    SimplePolygonMesh simpleMesh;
    Vector3 v_pos;
    vector<int> flags(0);
    for (Vertex v : mesh->vertices())
    {
        v_pos = geometry->inputVertexPositions[v];
        simpleMesh.vertexCoordinates.push_back(v_pos);
    }
    int id1;
    int id2;
    int id3;

    bool non_manifold = false;

    for (Face f : mesh->faces())
    {
        std::vector<size_t> polygon(3);

        Halfedge he = f.halfedge();
        id1 = he.vertex().getIndex();
        he = he.next();
        id2 = he.vertex().getIndex();
        he = he.next();
        id3 = he.vertex().getIndex();
        polygon[0] = id1;
        polygon[1] = id2;
        polygon[2] = id3;

        simpleMesh.polygons.push_back(polygon);
    }

    return simpleMesh;
}

std::vector<std::string> split(std::string s, std::string delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos)
    {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

Vector3 Get_bead_pos(std::string filename, int step)
{
    std::cout << "THe filename is" << filename << " \n";
    // Ok i need to get the vector for the position;
    std::ifstream Bead_data(filename);

    if (!Bead_data.is_open())
    {
        cerr << "Error opening the file! " << endl;
        return Vector3({0.0, 0.0, 0.0});
    }
    string line;
    int counter = 0;
    Vector3 Bead_pos({0.0, 0.0, 0.0});
    // the frame i want is (500*nsim)/100 = 5*nsim

    while (std::getline(Bead_data, line))
    {
        std::vector<std::string> splitted = split(line, " ");
        if (line[0] == '#')
        {
            continue;
        }

        if (counter == int(step / 50))
        {
            // We get the things we need
            // Bead_pos.x = std::stod(splitted[0]);
            Bead_pos = Vector3({std::stod(splitted[0]), std::stod(splitted[1]), std::stod(splitted[2])});
            break;
        }
        counter += 1;
    }
    if (counter < step / 50.0)
    {
        std::cout << "The output file doesnt get that far\n";
    }

    std::cout << "The initial bead position is " << Bead_pos.x << " " << Bead_pos.y << " " << Bead_pos.z << " \n";

    return Bead_pos;
}

void computeColors()
{

    Face_sizings = M3DG.Face_sizings();
    Vert_sizings = M3DG.Vert_sizing(Face_sizings);

    double F_sizing_max = 0.0;
    double V_sizing_max = 0.0;

    double F_sizing_min = 1e3;
    double V_sizing_min = 1e3;

    for (Vertex v : mesh->vertices())
    {
        if (Vert_sizings[v] > V_sizing_max)
            V_sizing_max = Vert_sizings[v];
        if (Vert_sizings[v] < V_sizing_min)
            V_sizing_min = Vert_sizings[v];
    }
    for (Face f : mesh->faces())
    {
        if (Face_sizings[f] > F_sizing_max)
            F_sizing_max = Face_sizings[f];
        if (Face_sizings[f] < F_sizing_min)
            F_sizing_min = Face_sizings[f];
    }

    // We have the min and the max
    for (Vertex v : mesh->vertices())
    {
        shaded_C.push_back({1.0, 0.45, 0.0});
    }
}

// Lets do our big one

VertexData<Vector3> EnergyGrad(std::string Energy, std::vector<double> energyConstants, int beadCount)
{
    // Can i just do a bunch of ifs
    if (Energy == "Bending")
        return Sim_handler.F_Bending(energyConstants);
    if (Energy == "Surface_tension")
        return Sim_handler.F_SurfaceTension(energyConstants);
    if (Energy == "Volume_constraint")
        return Sim_handler.F_Volume_constraint_2(energyConstants);
    if (Energy == "Area_constraint")
        return Sim_handler.F_Area_constraint(energyConstants);
    if (Energy == "Bending_tan")
        return Sim_handler.F_Bending_tan(energyConstants);
    if (Energy == "Edge_reg")
        return Sim_handler.F_Edge_reg(energyConstants);
    if (Energy == "Bead")
        return Sim_handler.Beads[beadCount]->Bead_I->Gradient();
    // If it is unkown we return a 0.0 vector
    VertexData<Vector3> Force(*mesh, Vector3{0.0, 0.0, 0.0});
    // For the bead energy i need to check how to do it
    return Force;
}

void updateTargetArea(size_t current_t, double origArea, double targetArea, double step, size_t req_steps)
{
    double currentTarget;
    if (req_steps < current_t)
    {
        currentTarget = targetArea;
    }
    else
    {
        // If not then
        currentTarget = origArea + step * current_t;
    }
    // Now i need to find it
    for (size_t i = 0; i < Energies.size(); i++)
    {
        if (Energies[i] == "Area_constraint")
        {

            Energy_constants[i][1] = currentTarget;
            Sim_handler.Energy_constants[i][1] = currentTarget;
            break;
        }
    }
}

double Integration_step(int timestep, bool Save, bool First_iter)
{
    double timeStep;
    M3DG.discreteTs = timestep;
    if (Save)
        Sim_data = std::ofstream(filename, std::ios_base::app);
    if (Integration == "Gradient_descent")
    {
        timeStep = M3DG.integrate(Sim_data, timestep, Bead_filenames, Save);
    }
    if (Integration == "BFGS")
    {
        M3DG.m = 10;
        timeStep = M3DG.integrate_BFGS(Sim_data, timestep, Bead_filenames, Save);
    }
    if (Integration == "BFGS-Normal")
    {
        M3DG.m = 10;
        if (First_iter)
        {
            std::cout << " Updating vertex normals for BFGS-Normal\n";
            M3DG.Sim_handler->update_vertex_normals();
            for (size_t i = 0; i < Beads.size(); i++)
            {
                Beads[i].state = "froze";
            }
        }
        // I will need the vertex normals (need something to know if they have been calculated)
        timeStep = M3DG.integrate_BFGS_Normal(Sim_data, timestep, Bead_filenames, Save);
    }
    if (Integration == "Newton")
    {
        // I need to know which are the constraints
        std::vector<std::string> Constraints;
        timeStep = M3DG.integrate_Newton(Sim_data, timestep, Bead_filenames, Save, Constraints, Constraints);
    }
    if (Integration == "Newton-Normal")
    {
        if (First_iter)
        {
            std::cout << "Beaut of contraints ";
            for (size_t i = 0; i < Constraints.size(); i++)
                std::cout << Constraints[i] << " ";
            for (size_t i = 0; i < Beads.size(); i++)
            {
                Beads[i].state = "froze";
            }
            Sim_handler.update_vertex_normals();
            Lagrange_mults.resize(Constraints.size());
            for (size_t i = 0; i < Constraints.size(); i++)
                Lagrange_mults(i) = 0.0f;
            Sim_handler.Lagrange_mult = Lagrange_mults;
            double A_target = 0.0;
            // Now i need to do the area constraint
            for (size_t i = 0; i < Energies.size(); i++)
            {
                if (Energies[i] == "Area_constraint")
                {
                    A_target = Energy_constants[i][1];
                    Sim_handler.Energy_constants[i][0] = 0.0;
                }
                if (Energies[i] == "Volume_constraint")
                {
                    Sim_handler.Energy_constants[i][0] = 0.0;
                    V_bar = Energy_constants[i][1];
                }
            }
            // Should i add the target volume and target area too?
            Sim_handler.Trgt_area = A_target;
            Sim_handler.Trgt_vol = V_bar;
            Sim_handler.Constraints = Constraints;
            // std::cout << "Solving the sytem to get the initial lagrange multipliers\n";
            // i WANT TO SOLVE THE rectangular system.
            Sim_handler.Calculate_Jacobian_Normal();
            Sim_handler.Calculate_gradient();
            // std::cout << "We called our gradients \n";
            Eigen::MatrixXd Jt = Sim_handler.Jacobian_constraints.transpose();
            Jt = Jt.block(0, 0, Sim_handler.mesh->nVertices(), 1);
            Eigen::VectorXd df(Sim_handler.mesh->nVertices());
            // std::cout << "The cuantities have been defined\n";
            std::cout << "The shape of Jt is" << Jt.rows() << " " << Jt.cols() << "\n";
            std::cout << "The shape of df is" << df.size() << "\n";
            for (size_t i = 0; i < Sim_handler.mesh->nVertices(); i++)
            {
                df(i) = dot(Sim_handler.Current_grad[i], Sim_handler.Vertex_normals[i]);
            }
            Eigen::MatrixXd J = Jt.transpose();
            // std::cout << "THe long awaited solve\n";
            Eigen::VectorXd delta_lambdas = (J * Jt).colPivHouseholderQr().solve(J * df);
            // std::cout << "Solved!\n";
            Sim_handler.Lagrange_mult(0) = delta_lambdas(0);
            // Sim_handler.Lagrange_mult(1) = delta_lambdas(1);
            // std::cout << "The initial lagrange multipliers are " << Sim_handler.Lagrange_mult.transpose() << "\n";
        }
        // Ok i need this thing to
        // I need to know if it is the first iteration or not ...
        timeStep = M3DG.integrate_Newton_Normal(Sim_data, timestep, Bead_filenames, Save, Constraints, Constraints);
        std::cout << "The timestep is " << timeStep << " \n";
    }
    if (Save)
        Sim_data.close();

    for (Bead b : Beads)
    {
        if (b.state == "manual")
        {
            b.update_state();
        }
    }

    return timeStep;
}

VertexData<Vector3> NewtonStep()
{

    for (size_t i = 0; i < Beads.size(); i++)
    {
        Beads[i].state = "froze";
    }
    std::cout << "Beaut of contraints ";
    for (size_t i = 0; i < Constraints.size(); i++)
        std::cout << Constraints[i] << " ";
    std::cout << "\n";
    Sim_handler.update_vertex_normals();
    Lagrange_mults.resize(Constraints.size());
    for (size_t i = 0; i < Constraints.size(); i++)
        Lagrange_mults(i) = 0.0f;
    Sim_handler.Lagrange_mult = Lagrange_mults;
    double A_target = 0.0;
    // Now i need to do the area constraint
    for (size_t i = 0; i < Energies.size(); i++)
    {
        if (Energies[i] == "Area_constraint")
        {
            A_target = Energy_constants[i][1];
            Sim_handler.Energy_constants[i][0] = 0.0;
        }
        if (Energies[i] == "Volume_constraint")
        {
            Sim_handler.Energy_constants[i][0] = 0.0;
            V_bar = Energy_constants[i][1];
        }
    }
    // Should i add the target volume and target area too?
    Sim_handler.Trgt_area = A_target;
    Sim_handler.Trgt_vol = V_bar;
    Sim_handler.Constraints = Constraints;
    // std::cout << "Solving the sytem to get the initial lagrange multipliers\n";
    // i WANT TO SOLVE THE rectangular system.
    Sim_handler.Calculate_Jacobian_Normal();
    Sim_handler.Calculate_gradient();
    // std::cout << "We called our gradients \n";
    Eigen::MatrixXd Jt = Sim_handler.Jacobian_constraints.transpose();
    Jt = Jt.block(0, 0, Sim_handler.mesh->nVertices(), 2);
    Eigen::VectorXd df(Sim_handler.mesh->nVertices());

    // std::cout << "The cuantities have been defined\n";
    // std::cout << "The shape of Jt is" << Jt.rows() << " " << Jt.cols() << "\n";
    // std::cout << "The shape of df is" << df.size() << "\n";
    for (size_t i = 0; i < Sim_handler.mesh->nVertices(); i++)
    {
        df(i) = dot(Sim_handler.Current_grad[i], Sim_handler.Vertex_normals[i]);
    }
    // std::cout << "THe long awaited solve\n";
    Eigen::VectorXd delta_lambdas = Jt.colPivHouseholderQr().solve(df);
    // std::cout << "Solved!\n";
    Sim_handler.Lagrange_mult(0) = delta_lambdas(0);
    Sim_handler.Lagrange_mult(1) = delta_lambdas(1);
    // std::cout << "The initial lagrange multipliers are " << Sim_handler.Lagrange_mult.transpose() << "\n";
    return M3DG.Newton_Normal_step(Sim_data, 0, Bead_filenames, false, Constraints, Constraints);
}

static const ImVec4 pressColor = ImColor::HSV(1. / 7.0f, 0.6f, 0.6f); // gold
static const ImVec4 releaseColor{0.35f, 0.61f, 0.49f, 0.62f};         // default green

static ImVec4 off_currColor = pressColor;
static ImVec4 V_sizing_currColor = releaseColor;
static ImVec4 F_sizing_currColor = releaseColor;

std::vector<ImVec4 *> sState = {&off_currColor, &V_sizing_currColor, &F_sizing_currColor};

int currFrame = 0;
bool autoPlaying = false;

void Callback_anim()
{
    frame_counter++;
    if (frame_counter >= 60)
        frame_counter = frame_counter % 60;
    ImGui::Text("The current frame is %i", frame_counter);

    bool updateFrameData = false;

    ImGui::Checkbox("Autoplay", &autoPlaying);
    if (autoPlaying && frame_counter % 5 == 0)
    {
        currFrame++;
        if (currFrame >= mesh_names.size())
        {
            currFrame = 0;
        }
        updateFrameData = true;
    }
    bool sliderUpdated = ImGui::SliderInt("Curr Frame", &currFrame, 0, mesh_names.size() - 1);
    if (currFrame >= mesh_names.size())
    {
        currFrame = 0;
    }
    updateFrameData = updateFrameData || sliderUpdated;
    ImGui::Text("Current mesh is %s", mesh_names[currFrame].c_str());
    if (updateFrameData && frame_counter % 10 == 0)
    {
        current_t = frames[currFrame];
        std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(basic_name + mesh_names[currFrame]);
        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();
        M3DG.mesh = mesh;
        M3DG.geometry = geometry;
        M3DG.Sim_handler->mesh = mesh;
        M3DG.Sim_handler->geometry = geometry;
        // Here
        for (size_t bi = 0; bi < Beads.size(); bi++)
        {
            Beads[bi].mesh = mesh;
            Beads[bi].geometry = geometry;
            Beads[bi].Bead_I->mesh = mesh;
            Beads[bi].Bead_I->geometry = geometry;
        }
        geometry->requireVertexPositions();
        geometry->inputVertexPositions = geometry->vertexPositions;
        psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
        M3DG.BFGS_iter = 0;
    }
    // Here go the display of the different gradients (easy to be said and to be done)
    if (ImGui::Button("Bending E grad"))
    {
        VertexData<Vector3> Bending_Egrad(*mesh, Vector3{0.0, 0.0, 0.0});
        Bending_Egrad = M3DG.Sim_handler->F_Bending(Energy_constants[0]);
        psMesh->addVertexVectorQuantity("Bending grad", Bending_Egrad);
    }
}

void Callback_qts()
{
    frame_counter++;
    if (frame_counter >= 60)
        frame_counter = frame_counter % 60;
    ImGui::Text("The current frame is %i", frame_counter);

    bool updateFrameData = false;
    bool currframeUpdated = false;
    bool integratorUpdate = false;

    if (integration_counter == integration_steps && integrating)
    {
        integrating = false;
        integratorUpdate = true;
        std::cout << "Done integrating, we will update the files available.\n";
        DIR *FD;
        struct dirent *in_file;
        const char *target_dir = basic_name.c_str();

        FD = opendir(target_dir);
        if (FD == NULL)
        {
            fprintf(stderr, "Error : Failed to open target directory - %s\n", strerror(errno));
            return;
        }
        frames.resize(0);
        /* Reading object (files, directories ...) from the folder) */
        while ((in_file = readdir(FD)))
        {

            if (DT_REG == in_file->d_type)
            {
                if (in_file->d_name[0] == 'm')
                {
                    // Ok SO THIS WILL GIVE ME
                    std::vector<std::string> split_result = split(in_file->d_name, "_");
                    // Ok now i want the second element
                    int frame = std::stoi(split_result[1]);
                    frames.push_back(frame);
                }
            }
        }
        std::sort(frames.begin(), frames.end());
        mesh_names.resize(0);
        for (int i = 0; i < frames.size(); i++)
        {
            mesh_names.push_back("membrane_" + std::to_string(frames[i]) + ".obj");
        }
        // We updated the size of the mesh names vector
        // Now i need to resize the BEAD data
        frameBeadPos.resize(0);
        std::vector<std::pair<size_t, Vector3>> BeadPosPair;
        std::ifstream BeadReadData;
        size_t frame;
        Vector3 beadPos;
        for (size_t bi = 0; bi < Beads.size(); bi++)
        {
            BeadPosPair.resize(0);
            BeadReadData = std::ifstream(basic_name + "Bead_" + std::to_string(bi) + "_data.txt");
            std::string line;
            while (std::getline(BeadReadData, line))
            {
                std::vector<std::string> splitted = split(line, " ");
                if (line[0] == '#')
                    continue;
                // If not i need to use it
                frame = std::stoi(splitted[0]);
                beadPos = Vector3{std::stod(splitted[1]), std::stod(splitted[2]), std::stod(splitted[3])};
                BeadPosPair.push_back(std::make_pair(frame, beadPos));
            }
            std::sort(BeadPosPair.begin(), BeadPosPair.end(), deterministic_sort);
            // now we have a list of frames that is being sorter
            frameBeadPos.push_back(BeadPosPair);
        }
        std::cout << "Filenames updated\n";
        currFrame = currFrame + savedSteps;
    }
    ImGui::Checkbox("Autoplay", &autoPlaying);
    if (autoPlaying && frame_counter % 10 == 0)
    {
        currFrame++;
        if (currFrame >= mesh_names.size())
        {
            currFrame = 0;
        }
        updateFrameData = true;
    }
    ImGui::SameLine(100);
    bool sliderUpdated = ImGui::SliderInt("Curr Frame", &currFrame, 0, mesh_names.size() - 1);
    currframeUpdated = currframeUpdated || ImGui::InputInt("Curremt Frame", &currFrame, 1, mesh_names.size() - 1);
    if (currFrame >= mesh_names.size())
    {
        currFrame = 0;
    }

    updateFrameData = updateFrameData || sliderUpdated;
    ImGui::Text("Current mesh is %s", mesh_names[currFrame].c_str());

    if ((updateFrameData && frame_counter % 5 == 0) || currframeUpdated || integratorUpdate)
    {
        // I need to load the mesh and then update it
        current_t = frames[currFrame];
        if (Area_constraint && updateArea)
            updateTargetArea(current_t, origA, A_bar, dA, minStepsReq);

        std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(basic_name + mesh_names[currFrame]);

        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();
        M3DG.mesh = mesh;
        M3DG.geometry = geometry;
        M3DG.Sim_handler->mesh = mesh;
        M3DG.Sim_handler->geometry = geometry;
        for (size_t bi = 0; bi < Beads.size(); bi++)
        {
            Beads[bi].mesh = mesh; // We reassing the mesh
            Beads[bi].geometry = geometry;
            Beads[bi].Bead_I->mesh = mesh;
            Beads[bi].Bead_I->geometry = geometry;
        }
        geometry->requireVertexPositions();

        psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
        polyscope::state::facePickIndStart = mesh->nVertices();
        polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
        polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

        if (!integratorUpdate)
        {
            last_remesh = 0;
            M3DG.BFGS_iter = 0;
            continueLastCall = false;
        }

        CurrentBeadPos.resize(0);

        for (size_t bi = 0; bi < Beads.size(); bi++)
        {
            CurrentBeadPos.push_back(frameBeadPos[bi][currFrame].second);
            Beads[bi].Pos = CurrentBeadPos[bi];
        }
        if (Beads.size() > 0)
        {
            psCloud = polyscope::registerPointCloud("Great beads", CurrentBeadPos);
            psCloud->setPointRadius(Beads[0].sigma, false); // TODO adjust o different sizes
        }
    }

    ImGuiTreeNodeFlags flag2 = ImGuiTreeNodeFlags_None;

    ImGuiTreeNodeFlags flag = ImGuiTreeNodeFlags_None;
    if (ImGui::TreeNodeEx("Integration menu", flag))
    {
        // Ok so i am here and I already decided there are two things i need to do
        // One is the integration settings and one is the actual integration
        ImGui::Text(("The current integration method is " + Integration).c_str());
        if (ImGui::Button("Integrate!"))
        {
            integrating = true;
            savedSteps = 0;
            integration_counter = 0;
            // last_remesh = 0;
            current_t = frames[currFrame];
            std::cout << "We are starting from t = " << current_t << " \n";
        }
        ImGui::SameLine(150);
        ImGui::Checkbox("Continue last call", &continueLastCall);

        if (ImGui::TreeNodeEx("Options"))
        {
            // int integrationSteps = integtration_steps;
            // int remeshEvery;

            ImGui::InputInt("Integration steps", &integration_steps, 10, 1000);
            ImGui::InputInt("Remesh every", &remesh_every, 1, 50);
            ImGui::InputInt("Save every", &save_interval);
            ImGui::Checkbox("Update target area", &updateArea);
            const char *integrationMethods[] = {"Gradient_descent", "BFGS", "BFGS-Normal", "Newton", "Newton-Normal"};
            const string integrationMethodsString[] = {"Gradient_descent", "BFGS", "BFGS-Normal", "Newton", "Newton-Normal"};
            static int currentIntegration = 1;
            const char *preview = integrationMethods[currentIntegration];
            if (ImGui::BeginCombo("Integration Method", preview))
            {
                for (int n = 0; n < IM_ARRAYSIZE(integrationMethods); n++)
                {
                    bool isSelected = (currentIntegration == n);
                    if (ImGui::Selectable(integrationMethods[n], isSelected))
                    {
                        currentIntegration = n;
                        Integration = integrationMethodsString[n];
                        // Apply weapon change to game
                    }
                    if (isSelected)
                    {
                        ImGui::SetItemDefaultFocus();
                        //
                    }
                }
                ImGui::EndCombo();
            }
            ImGui::TreePop();
        }
        if (ImGui::TreeNodeEx("Newton funoptions"))
        {
            // 3 THINGS
            // 1. Calculate the newton step
            // 2. stepsize
            // 3. step that amount

            if (ImGui::Button("Calculate Newton step"))
            {
                NewtonStepSaved = NewtonStep();
                psMesh->addVertexVectorQuantity("Newton step", NewtonStepSaved);
            }

            ImGui::InputFloat("Step size", &NewtonTs, 1.0, 1.0, 4);
            if (ImGui::Button("Take Newton step"))
            {
                geometry->vertexPositions += NewtonTs * NewtonStepSaved;
                psMesh->updateVertexPositions(geometry->vertexPositions);
            }
            ImGui::TreePop();
        }

        ImGui::TreePop();
    }
    // How often do we check the small angle? every 10?

    if (integrating && frame_counter % 2 == 0)
    {

        // Ok now i can integrate every second frame(hopefully this helps us)
        integration_counter += 1;
        bool Save = integration_counter % save_interval == 0;
        // I should save before running
        std::cout << "We are at t = " << current_t + integration_counter << " \n";

        // i should add the small angle flag
        bool flagSmallAngle = false;
        if (integration_counter % 5 == 0)
        {
            geometry->requireCornerAngles();
            for (Corner c : mesh->corners())
            {
                if (geometry->cornerAngles[c] < 0.2)
                {
                    flagSmallAngle = true;
                    break;
                }
            }
            if (flagSmallAngle)
            {
                remeshSmallAngles(*mesh, *geometry, Options);
                M3DG.BFGS_iter = 0;
                mesh->compress();
                geometry->refreshQuantities(); // i guess yes
            }
            geometry->unrequireCornerAngles();
            // we do the small angle check
        }
        if (((integration_counter - last_remesh) >= remesh_every && remesh_every > 0) || flagSmallAngle || M3DG.remesh_flag)
        {
            M3DG.remesh_flag = false;
            last_remesh = integration_counter;
            int n_op = remesh(*mesh, *geometry, Options);
            // Here i need to update the remesh_every hehe
            if (adapt_remesh)
            {
                // HERE WE DO THE THINGY
                if (n_op > 50)
                {
                    remesh_every = remesh_every / 2;
                    if (remesh_every < 1)
                    {
                        remesh_every = 1;
                    }
                }
                else if (n_op < 20)
                {
                    remesh_every = remesh_every + 10;
                }
                else
                {
                    remesh_every = remesh_every + 1;
                }
                if (remesh_every > 200)
                {
                    remesh_every = 200;
                }
            }
            if (Count_remesh)
            {
                Remeshing_count = std::ofstream(basic_name + "Remeshing_count.txt", std::ios_base::app);
                Remeshing_count << current_t << " " << n_op << " " << remesh_every << " " << -1 << "\n";
                Remeshing_count.close();
            }
            mesh->compress();
            geometry->refreshQuantities();
            M3DG.BFGS_iter = 0;
            // I need to remesh
        }
        if (Area_constraint && updateArea)
            updateTargetArea(current_t + integration_counter, origA, A_bar, dA, minStepsReq);
        double timestepIter = Integration_step(integration_counter + current_t, Save, integration_counter == 1 && !continueLastCall);

        if (Save)
        {
            savedSteps += 1;
            Save_mesh(basic_name, current_t + integration_counter);
        }
        if (timestepIter <= 0)
        {
            std::cout << "Sim is not progresssing we will stop integrating\n";
            integration_counter = integration_steps;
        }
    }

    if (ImGui::TreeNodeEx("Remeshing menu", flag))
    {
        if (ImGui::IsItemClicked())
        {
        }
        if (ImGui::SmallButton("remesh"))
        {
            // I WANT TO REMESH
            polyscope::state::subset.vertices.clear();
            polyscope::state::subset.edges.clear();
            polyscope::state::subset.faces.clear();
            int n_op;
            n_op = remesh(*mesh, *geometry, Options);
            mesh->compress();
            geometry->refreshQuantities();
            // M3DG.mesh = mesh;
            // M3DG.geometry = geometry;
            psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
            polyscope::state::facePickIndStart = mesh->nVertices();
            polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
            polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();
        }
        ImGui::SameLine(100);
        if (ImGui::SmallButton("Smoothing"))
        {
            double move = smoothByLaplacian(*mesh, *geometry);
            std::cout << "The movement is " << move << " \n";
            psMesh->updateVertexPositions(geometry->vertexPositions);
        }
        // I want to add a button for triangle quality.
        if (ImGui::Button("^Qual"))
        {
            geometry->requireCornerAngles();
            FaceData<double> min_angle(*mesh, 0.0);
            double minAngle;
            for (Face f : mesh->faces())
            {
                //
                minAngle = 3.0;
                for (Corner c : f.adjacentCorners())
                {
                    if (geometry->cornerAngles[c] < minAngle)
                        minAngle = geometry->cornerAngles[c];
                }
                if (minAngle > angleThresh)
                    minAngle = 1.0;
                min_angle[f] = minAngle;
            }
            psMesh->addFaceScalarQuantity("Min angle", min_angle);
            geometry->unrequireCornerAngles();
        }
        ImGui::SameLine(100);
        if (ImGui::Button("Remesh Small"))
        {
            Options.angleThresh = angleThresh;
            remeshSmallAngles(*mesh, *geometry, Options);
            geometry->refreshQuantities();
            psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
            polyscope::state::facePickIndStart = mesh->nVertices();
            polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
            polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();
        }
        ImGui::SameLine(200);
        ImGui::DragFloat("Min Angle", &angleThresh, 0.05f, 0.05f, 0.5f, "%.2f");

        if (ImGui::TreeNodeEx("RemeshingOptions"))
        {
            if (ImGui::Button("Vertex Sizings"))
            {
                // geometry->refreshQuantities();
                geometry->purgeQuantities();
                double refine_angle = Options.refine_angle;
                double min_absolute_length = Options.min_absolute_length;
                double max_absolute_length = Options.max_absolute_length;
                geometry->requireFaceSizing();
                for (Face f : mesh->faces())
                {
                    geometry->faceSizing[f] = clamp(geometry->faceSizing[f] / (refine_angle * refine_angle),
                                                    1.0 / (max_absolute_length * max_absolute_length),
                                                    1.0 / (min_absolute_length * min_absolute_length));
                }
                psMesh->addFaceScalarQuantity("Face_sizings", geometry->faceSizing);
                std::cout << "Requiring the evertex sizing\n";
                geometry->requireVertexSizing();

                // std::cout << "THe first vertex sizing is after " << Vert_sizings[0] << "\n";
                psMesh->addVertexScalarQuantity("Vertex sizings", geometry->vertexSizing);
                // geometry->refreshQuantities();
                geometry->unrequireFaceSizing();
                geometry->unrequireVertexSizing();
            }
            if (ImGui::IsItemClicked())
            {

                // Mark rendered node as being clicked
            }

            static float sizeMin = Options.min_absolute_length;
            static float sizeMax = Options.max_absolute_length;
            static float aspectMin = Options.aspect_min;
            static float refineAngle = Options.refine_angle;
            ImGui::DragFloat("Min size", &sizeMin, 0.001f, 0.001f, 0.1f, "%.3f");
            Options.min_absolute_length = sizeMin;
            ImGui::DragFloat("Max size", &sizeMax, 0.05f, 0.05f, 1.0f, "%.2f");
            Options.max_absolute_length = sizeMax;
            ImGui::DragFloat("Aspect min", &aspectMin, 0.1f, 0.01f, 0.5f, "%.2f");
            Options.aspect_min = aspectMin;
            ImGui::DragFloat("Refine Angle", &refineAngle, 0.3f, 1.0f, 0.1f, "%.1f");
            Options.refine_angle = refineAngle;
            // This is where we will add the different options to change
            ImGui::TreePop();
        }
        if (ImGui::Button("Display normals"))
        {
            M3DG.Sim_handler->update_vertex_normals();
            psMesh->addVertexVectorQuantity("Normals", M3DG.Sim_handler->Vertex_normals);
        }

        ImGui::TreePop(); // This is required at the end of the if block
    }

    ImGuiTreeNodeFlags flag3 = ImGuiTreeNodeFlags_None;
    // Ok i want to do the energies
    // so the thing is that every energy has a different function that gets called. And i want the correct one to be called.
    if (ImGui::TreeNodeEx("Energy menu", flag))
    {
        // Here we show the full gradient too jeje
        if (ImGui::Button("Display gradient"))
        {
            Sim_handler.Calculate_gradient();
            psMesh->addVertexVectorQuantity("Full gradient", Sim_handler.Current_grad);
        }
        ImGui::SameLine(150);
        if (ImGui::Button("Display Newton step"))
        {
            psMesh->addVertexVectorQuantity("Newton step", NewtonStep());
        }

        static std::vector<std::vector<float>> Menu_E_constants(Og_Energy_constants.size());
        // Ok so now
        for (size_t i = 0; i < Energy_constants.size(); i++)
        {
            std::vector<float> values;
            for (size_t j = 0; j < Energy_constants[i].size(); j++)
            {
                values.push_back(float(Energy_constants[i][j]));
            }
            Menu_E_constants[i] = values;
        }

        if (ImGui::IsItemClicked())
        {

            // Mark rendered node as being clicked
        }
        // Ok so now we do every energy as a treenode
        int beadCount = 0;
        for (size_t i = 0; i < Energies.size(); i++)
        {
            std::string E_name = Energies[i];
            if (E_name == "Bead")
                continue;
            if (ImGui::TreeNodeEx(E_name.c_str()))
            {
                // Here we add the button to show the gradient
                if (ImGui::Button("Show Gradient"))
                {
                    psMesh->addVertexVectorQuantity(E_name + " gradient", EnergyGrad(E_name, Energy_constants[i], beadCount));
                }
                for (size_t j = 0; j < Menu_E_constants[i].size(); j++)
                {
                    // static float Value = Energy_constants[i][0];
                    ImGui::DragFloat(("Constant " + std::to_string(j)).c_str(), &Menu_E_constants[i][j], Og_Energy_constants[i][j] * 0.1, Og_Energy_constants[i][j] * 10.0, Og_Energy_constants[i][j] * 0.1, "%.3f");
                    Energy_constants[i][j] = Menu_E_constants[i][j];
                    // Sim_handler.Energy_constants[i][j] = Energy_constants[i][j];
                }

                ImGui::TreePop();
            }
            // OK now i will add just a button to show the gradient of the beads
            // In theory i need something that changes the parameters but thats too much man
        }
        for (size_t bi = 0; bi < Beads.size(); bi++)
        {
            // we want to add the energy quantity here
            if (ImGui::TreeNodeEx(("Bead " + std::to_string(bi)).c_str()))
            {
                if (ImGui::Button("Show Gradient"))
                {
                    psMesh->addVertexVectorQuantity("Bead " + std::to_string(bi) + " gradient", EnergyGrad("Bead", Energy_constants[0], bi));
                }

                ImGui::TreePop();
            }
        }

        if (Beads.size() > 0)
        {
            // In case there are beads
            if (ImGui::Button("Show bead forces"))
            {
                // IN this case i will add the bead forces
                CurrentBeadForce.resize(0);
                for (size_t bi = 0; bi < Beads.size(); bi++)
                {
                    std::cout << "The bead is located at " << Beads[bi].Pos << " \n";
                    Beads[bi].Bead_I->mesh = mesh;
                    Beads[bi].Bead_I->geometry = geometry;

                    Beads[bi].Bead_I->Gradient();
                    CurrentBeadForce.push_back(Beads[bi].Total_force);
                    std::cout << "The bead force at 0 is " << CurrentBeadForce[bi] << " \n";
                }
                psCloud->addVectorQuantity("Bead Forces", CurrentBeadForce);
            }
        }
        ImGui::TreePop();
    }
    if (ImGui::TreeNodeEx("Transformation menu", flag))
    {
        if (ImGui::TreeNodeEx("Rescale"))
        {
            // THere are 3 factors Sx Sy and Sz
            ImGui::InputFloat("Scale x", &Sx, 0.1, 1.0, 4);
            ImGui::InputFloat("Scale y", &Sy, 0.1, 1.0, 4);
            ImGui::InputFloat("Scale z", &Sz, 0.1, 1.0, 4);

            if (ImGui::Button("Rescale!"))
            {
                Vector3 scaleVector = Vector3{Sx, Sy, Sz};
                geometry->rescale(scaleVector);
                psMesh->updateVertexPositions(geometry->inputVertexPositions);
                geometry->vertexPositions = geometry->inputVertexPositions;
                polyscope::requestRedraw();
            }
            // I need to refresh no

            ImGui::TreePop();
        }

        ImGui::TreePop();
    }

    if (ImGui::TreeNodeEx("Selection menu", flag))
    {
        if (ImGui::Button("Show selection"))
        {
            showSelected();
        }
        ImGui::SameLine(150);
        if (ImGui::Button("Reset selection"))
        {
            polyscope::state::subset.vertices.clear();
            polyscope::state::subset.edges.clear();
            polyscope::state::subset.faces.clear();
            showSelected();
            // redraw();
        }
        ImGui::SameLine(280);
        if (ImGui::Button("Remesh_selection"))
        {
            // I want to get all the faces
            Face f;
            std::vector<Face> active_faces;
            for (std::set<size_t>::iterator it = polyscope::state::subset.faces.begin();
                 it != polyscope::state::subset.faces.end(); ++it)
            {
                f = mesh->face(*it);
                active_faces.push_back(f);
            }
            remesh(*mesh, *geometry, Options, active_faces);
            geometry->requireVertexPositions();
            geometry->requireVertexIndices();
            geometry->requireEdgeIndices();

            mesh->compress();
            M3DG.mesh = mesh;
            M3DG.geometry = geometry;
            psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());

            polyscope::state::facePickIndStart = mesh->nVertices();
            polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
            polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();
        }

        ImGui::TreePop();
    }
}

void Callback_curvatures()
{
    double HB = Sim_handler.E_Bending(Energy_constants[0]);
    double HL = Sim_handler.E_Laplace(Energy_constants[0]);
    double EC = 4 * 3.1415926535 * Energy_constants[0][0];
    ImGui::Text("The Bending Energy is %0.1f ", HB);
    ImGui::Text("The Laplace Energy is %0.1f ", HL);
    ImGui::Text("The expected Energy is %0.1f", EC);

    if (ImGui::Button("Mean curvature"))
    {
        VertexData<double> Mean_curv(*mesh, 0.0);
        for (Vertex v : mesh->vertices())
        {
            Mean_curv[v] = geometry->scalarMeanCurvature(v) / geometry->barycentricDualArea(v);
        }
        psMesh->addVertexScalarQuantity("Mean curvature dihedrals", Mean_curv);
    }
    if (ImGui::Button("Mean curvature 2"))
    {
        VertexData<double> Mean_curv_2(*mesh, 0.0);
        for (Vertex v : mesh->vertices())
        {
            Mean_curv_2[v] = geometry->vertexNormalMeanCurvature(v).norm() / geometry->barycentricDualArea(v);
        }
        psMesh->addVertexScalarQuantity("Norm of normal cotangen", Mean_curv_2);
    }
    if (ImGui::Button("Bending E grad"))
    {
        VertexData<Vector3> Bending_Egrad(*mesh, Vector3{0.0, 0.0, 0.0});
        Bending_Egrad = M3DG.Sim_handler->F_Bending(Energy_constants[0]);
        psMesh->addVertexVectorQuantity("Bending Egrad", Bending_Egrad);
    }
    if (ImGui::Button("Laplace E grad"))
    {
        VertexData<Vector3> Laplace_Egrad(*mesh, Vector3{0.0, 0.0, 0.0});
        Laplace_Egrad = M3DG.Sim_handler->F_Laplace(Energy_constants[0]);
        psMesh->addVertexVectorQuantity("Laplace Egrad", Laplace_Egrad);
    }
}

void functionCallback()
{

    // std::cout<<"Functioncallback\n";
    ImGui::Text("Total angle defect: %0.1fpi", TOTAL_ANGLE_DEFECT / M_PI);
    ImGui::Text("Euler characteristic: %zu", EULER_CHARACTERISTIC);
    ImGui::Text("Saved slots %i", Save_slot - 1);

    for (size_t i = 0; i < Beads.size(); i++)
    {
        ImGui::Text("The bead position is %0.2f %0.2f %0.2f", Beads[i].Pos.x, Beads[i].Pos.y, Beads[i].Pos.z);
    }
    ImGui::Text("");

    ImGui::PushStyleColor(ImGuiCol_Button, *sState[0]);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, *sState[0]);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, *sState[0]);

    // if(Energies[0]=="Laplace"){
    if (ImGui::Button("Laplace"))
    {
        Gradient_vertex = Sim_handler.F_Laplace(Energy_constants[0]);

        std::cout << "Calculating Laplace\n";
        psMesh->addVertexVectorQuantity("Laplace", Gradient_vertex);
    }
    if (ImGui::Button("Edge_reg"))
    {
        Gradient_vertex = Sim_handler.F_Edge_reg(Energy_constants[2]);

        std::cout << "Calculating Edge reg energy\n";
        psMesh->addVertexVectorQuantity("Edge reg", Gradient_vertex);
    }
    if (ImGui::Button("Edge_reg_2"))
    {
        Gradient_vertex = Sim_handler.F_Edge_reg_2(Energy_constants[0]);
        std::cout << "Calculating Edge reg energy 2\n";
        psMesh->addVertexVectorQuantity("Edge reg 2", Gradient_vertex);
    }
    if (ImGui::Button("Bending"))
    {
        Gradient_vertex = Sim_handler.F_Bending(Energy_constants[0]);
        std::cout << "Calculating Bending\n";
        psMesh->addVertexVectorQuantity("Bending", Gradient_vertex);
    }
    if (ImGui::Button("Surface tension"))
    {
        Gradient_vertex = Sim_handler.F_SurfaceTension(Energy_constants[1]);
        std::cout << "Calculating Surface tension\n";
        psMesh->addVertexVectorQuantity("Surface tension", Gradient_vertex);
    }
    if (ImGui::Button("Bead"))
    {
        Gradient_vertex = Sim_handler.Beads[0]->Bead_I->Gradient();
        std::cout << "Calculating Bead\n";
        psMesh->addVertexVectorQuantity("Bead", Gradient_vertex);
    }
    if (ImGui::Button("Volume"))
    {
        Gradient_vertex = Sim_handler.F_Volume_constraint_2(Energy_constants[3]);
        std::cout << "Calculating Volume constraint\n";
        psMesh->addVertexVectorQuantity("Volume constraint", Gradient_vertex);
    }
    // }

    if (ImGui::Button("Display grad"))
    {
        // Results/Debug_remesh_trial/Bending_1.0000_Bead_radius_0.3000_str_400.0000_Nsim_8
        std::cout << "Calculating gradient\n";
        M3DG.Sim_handler->Calculate_gradient();
        Gradient_vertex = M3DG.Sim_handler->Current_grad;
        // Gradient_vertex = M3DG.Bending(0.0) + Beads[0].Gradient();
        psMesh->addVertexVectorQuantity("Gradient", Gradient_vertex);
    }

    if (Eigenvectors_Bending.size() > 0)
    {

        if (ImGui::Button("First eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[0][3 * v.getIndex()],
                                         Eigenvectors_Bending[0][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[0][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(0), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 1)
    {

        if (ImGui::Button("Second eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[1][3 * v.getIndex()],
                                         Eigenvectors_Bending[1][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[1][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(1), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 2)
    {

        if (ImGui::Button("Third eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[2][3 * v.getIndex()],
                                         Eigenvectors_Bending[2][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[2][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(2), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 3)
    {

        if (ImGui::Button("Forth eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[3][3 * v.getIndex()],
                                         Eigenvectors_Bending[3][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[3][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(3), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 4)
    {

        if (ImGui::Button("Fifth eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[4][3 * v.getIndex()],
                                         Eigenvectors_Bending[4][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[4][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(4), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 5)
    {

        if (ImGui::Button("Sixth eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[5][3 * v.getIndex()],
                                         Eigenvectors_Bending[5][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[5][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(5), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 6)
    {

        if (ImGui::Button("Seventh eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[6][3 * v.getIndex()],
                                         Eigenvectors_Bending[6][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[6][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(6), Eigenvector);
        }
    }
    if (Eigenvectors_Bending.size() > 7)
    {

        if (ImGui::Button("Eigth eigenval "))
        {
            VertexData<Vector3> Eigenvector(*mesh, Vector3{0.0, 0.0, 0.0});
            for (Vertex v : mesh->vertices())
            {
                Eigenvector[v] = Vector3{Eigenvectors_Bending[7][3 * v.getIndex()],
                                         Eigenvectors_Bending[7][3 * v.getIndex() + 1],
                                         Eigenvectors_Bending[7][3 * v.getIndex() + 2]};
            }
            std::cout << "Displaying bending eigenvals\n";
            psMesh->addVertexVectorQuantity("Bending eigenvals " + std::to_string(7), Eigenvector);
        }
    }

    if (ImGui::InputInt("Timestep to inspect", &inspect_timestep))
        ;

    if (ImGui::Button("Display grad to inspect"))
    {

        VertexData<Vector3> Gradient;
        VertexData<Vector3> Gradient2;
        VertexData<Vector3> Gradient3;
        Vector3 Sample;
        std::cout << "Inspecting timestep " << inspect_timestep << "\n";

        std::string filename = basic_name + "RHS_Norm" + std::to_string(inspect_timestep) + ".txt";

        std::ifstream f(filename);
        string s;
        if (!f.is_open())
        {
            std::cout << "The file " << filename << " does not exist\n";
        }
        else
        {
            getline(f, s);
            std::vector<std::string> splitted = split(s, " ");
            if (splitted.size() > 100)
            {
                std::cout << "We are loading the memebrane now \n";
                std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(basic_name + "membrane_" + std::to_string(inspect_timestep) + ".obj");
                mesh = mesh_uptr.release();
                geometry = geometry_uptr.release();
                M3DG.mesh = mesh;
                M3DG.geometry = geometry;
                Sim_handler.mesh = mesh;
                Sim_handler.geometry = geometry;

                std::cout << "The mesh has " << mesh->nVertices() << " vertices and " << mesh->nFaces() << " faces\n";
                geometry->requireVertexPositions();
                std::cout << "Registering surfaces\n";
                Gradient_vertex = VertexData<Vector3>(*mesh, Vector3{0.0, 0.0, 0.0});

                Gradient = VertexData<Vector3>(*mesh, Vector3{0.0, 0.0, 0.0});
                Gradient2 = VertexData<Vector3>(*mesh, Vector3{0.0, 0.0, 0.0});
                Gradient3 = VertexData<Vector3>(*mesh, Vector3{0.0, 0.0, 0.0});

                // We are going to create the file
                // std::cout<<"The line is "<< s <<"\n";

                std::cout << "The size of splitted is " << splitted.size() << "\n";
                std::cout << "THe number of vertices times 3 is " << mesh->nVertices() * 3 << "\n";
                for (size_t i = 0; i < mesh->nVertices(); i++)
                {
                    // std::cout<<"The splitted size is "<< splitted.size() <<"\n"
                    // std::cout<<"The next line is " << splitted[3*i] << " " << splitted[3*i+1] << " " << splitted[3*i+2] <<"\n";
                    // std::cout<<"The size of splitted is "<< splitted.size() <<"\n";
                    // std::cout<<"The number of vertices is \n";
                    Sample.x = std::stod(splitted[3 * i]);
                    Sample.y = std::stod(splitted[3 * i + 1]);
                    Sample.z = std::stod(splitted[3 * i + 2]);
                    Gradient[i] = Sample;
                }
            }
            getline(f, s);
            splitted = split(s, " ");
            for (size_t i = 0; i < mesh->nVertices(); i++)
            {
                // std::cout<<"The next line is " << splitted[3*i] << " " << splitted[3*i+1] << " " << splitted[3*i+2] <<"\n";
                Sample.x = std::stod(splitted[3 * i]);
                Sample.y = std::stod(splitted[3 * i + 1]);
                Sample.z = std::stod(splitted[3 * i + 2]);
                Gradient2[i] = Sample;
            }
            getline(f, s);
            splitted = split(s, " ");
            for (size_t i = 0; i < mesh->nVertices(); i++)
            {
                // std::cout<<"The next line is " << splitted[3*i] << " " << splitted[3*i+1] << " " << splitted[3*i+2] <<"\n";
                Sample.x = std::stod(splitted[3 * i]);
                Sample.y = std::stod(splitted[3 * i + 1]);
                Sample.z = std::stod(splitted[3 * i + 2]);
                Gradient3[i] = Sample;
            }

            psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
            psMesh->addVertexVectorQuantity("Gradient Lag", Gradient);
            psMesh->addVertexVectorQuantity("Gradient E", Gradient2);
            psMesh->addVertexVectorQuantity("Jacobian ", Gradient3);
            // Now that we have the gradient we need to display it
        }
    }
    int bead_counter = 0;

    if (ImGui::Button("Display newton"))
    {
        std::cout << "Calculating gradient newton\n";
        VertexData<Vector3> Gradient_newton;
        if (first_newton)
        {
            Eigen::VectorXd Lagrange_mults(0);

            Sim_handler.Lagrange_mult = Lagrange_mults;
            // Sim_handler.Constraints = std::vector<std::string>{"Volume"};
            Sim_handler.Constraints = std::vector<std::string>{};
            first_newton = false;
        }
        // M3DG.integrate_Newton(Sim_data,0.0,Energies,false,std::vector<std::string>{"Volume","CMx","CMy","CMz","Rx","Ry","Rz"},std::vector<std::string>{"Files"});
        M3DG.integrate_Newton(Sim_data, 0.0, Energies, false, std::vector<std::string>{}, std::vector<std::string>{"Files"});

        Gradient_newton = Sim_handler.Current_grad;

        psMesh->addVertexVectorQuantity("Gradient Newton", Gradient_newton);
        redraw();
    }
    if (ImGui::Button("Step newton"))
    {
        std::cout << "Calculating gradient newton\n";
        VertexData<Vector3> Gradient_newton;
        if (first_newton)
        {
            Eigen::VectorXd Lagrange_mults(0);
            Sim_handler.Lagrange_mult = Lagrange_mults;
            Sim_handler.Constraints = std::vector<std::string>{};
            // Sim_handler.Constraints = std::vector<std::string>{"Volume","Area","CMx","CMy","CMz","Rx","Ry","Rz"};
            first_newton = false;
        }

        // M3DG.integrate_Newton(Sim_data,0.0,Energies,false,std::vector<std::string>{"Volume","Area","CMx","CMy","CMz","Rx","Ry","Rz"},std::vector<std::string>{"Files"});
        M3DG.integrate_Newton(Sim_data, 0.0, Energies, false, std::vector<std::string>{}, std::vector<std::string>{"Files"});

        // Gradient_newton = Sim_handler.Current_grad;
        // geometry->inputVertexPositions += Gradient_newton;
        std::cout << "Redrawing\n";
        redraw();
        // psMesh->addVertexVectorQuantity("Gradient Newton", Gradient_newton);
    }

    // if(ImGui::Button(Energies[1].c_str()))
    // {
    //     VertexData<Vector3> Gradient_temp = M3DG.Sim_handler->F_Bending(Energy_constants[1]);
    //     psMesh->addVertexVectorQuantity("Gradient Bending", Gradient_temp);

    // }

    // if(ImGui::Button(Energies[0].c_str()))
    // {
    //     VertexData<Vector3> Gradient_temp = M3DG.Sim_handler->F_SurfaceTension(Energy_constants[0]);
    //     psMesh->addVertexVectorQuantity("Gradient Surface tension", Gradient_temp);
    // }

    // if(Beads.size()>0){

    // // VertexData<Vector3> Gradient_temp = Beads[bead_counter].Gradient();
    // if(ImGui::Button(Energies[2].c_str()))
    // {
    //     VertexData<Vector3> Gradient_temp = Beads[0].Gradient();
    //     std::cout<<"The gradient of this bead is " << Gradient_temp[0].x << " " << Gradient_temp[0].y << " " << Gradient_temp[0].z << "\n";
    //     psMesh->addVertexVectorQuantity("Gradient Bead 1" , Gradient_temp);

    // }
    // if(ImGui::Button("Bead 2"))
    // {
    //     VertexData<Vector3> Gradient_temp = Beads[1].Gradient();
    //     std::cout<<"The gradient of this bead is " << Gradient_temp[0].x << " " << Gradient_temp[0].y << " " << Gradient_temp[0].z << "\n";
    //     psMesh->addVertexVectorQuantity("Gradient Bead 2" , Gradient_temp);
    // }
    // }

    if (ImGui::Button("Save current state"))
    {

        Saved_mesh = save_geometry(mesh, geometry);
        Save_mesh(basic_name, Save_slot);
        if (Beads.size() > 0)
        {

            std::ofstream bead_saved(basic_name + "Bead_data_0.txt", ios_base::app);
            bead_saved << Beads[0].Pos.x << " " << Beads[0].Pos.y << " " << Beads[0].Pos.z << "\n";
            bead_saved.close();
            Saved_beadpos[0] = Beads[0].Pos;
        }

        Save_slot++;
        std::cout << "SAVE SLOT IS " << Save_slot << " \n";
    }
    if (ImGui::Button("Reload saved state"))
    {

        std::cout << "reloading\n";
        delete mesh;
        delete geometry;
        std::cout << "deleted\n";

        auto lvals = makeManifoldSurfaceMeshAndGeometry(Saved_mesh.polygons, Saved_mesh.vertexCoordinates);

        std::tie(mesh_uptr, geometry_uptr) = std::tuple<std::unique_ptr<ManifoldSurfaceMesh>,
                                                        std::unique_ptr<VertexPositionGeometry>>(std::move(std::get<0>(lvals)), // mesh
                                                                                                 std::move(std::get<1>(lvals)));

        std::cout << "reassigned\n";

        mesh = mesh_uptr.release();
        geometry = geometry_uptr.release();

        M3DG.mesh = mesh;
        M3DG.geometry = geometry;
        std::cout << "deleted\n";
        if (Beads.size() > 0)
        {
            Beads[0].Pos = Saved_beadpos[0];
        }
        geometry->requireVertexPositions();
        std::cout << "relocated\n";
        psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
        std::cout << "re registered\n";
    }
    ImGui::InputInt("Load slot", &Load_slot);
    if (ImGui::Button("Reload saved slot"))
    {
        std::cout << "1\n";
        std::string filepath = basic_name + "membrane_" + std::to_string(Load_slot) + ".obj";
        std::cout << "2\n";
        delete mesh;
        std::cout << "3\n";
        delete geometry;
        std::cout << "4\n";
        std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
        std::cout << "5\n";

        std::ifstream bdata(basic_name + "Bead_data_0.txt");
        int counter = 0;
        string s;
        Vector3 NewPos;
        while (getline(bdata, s))
        {
            //
            if (counter == Load_slot)
            {
                std::vector<std::string> splitted_line = split(s, " ");

                NewPos.x = std::stod(splitted_line[0]);
                NewPos.y = std::stod(splitted_line[1]);
                NewPos.z = std::stod(splitted_line[2]);
                break;
            }
            counter += 1;
        }
        std::cout << "New pos is " << NewPos << " \n";

        mesh = mesh_uptr.release();
        std::cout << "6\n";
        geometry = geometry_uptr.release();
        std::cout << "7\n";
        M3DG.mesh = mesh;
        std::cout << "8\n";

        M3DG.geometry = geometry;
        std::cout << "9\n";
        geometry->requireVertexPositions();
        Beads[0].Pos = NewPos;

        std::cout << "10\n";
        psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());

        redraw();
    }

    ImGui::InputDouble("Scaling_factor", &scaling_factor);
    if (ImGui::Button("Rescale mesh "))
    {
        for (Vertex v : mesh->vertices())
        {
            geometry->inputVertexPositions[v] *= scaling_factor;
        }

        geometry->refreshQuantities();
        psMesh->updateVertexPositions(geometry->inputVertexPositions);
    }

    if (ImGui::Button("Step flow"))
    {
        for (size_t i = 0; i < integration_steps; i++)
        {

            double dt = M3DG.integrate(Sim_data, 0.0, Bead_filenames, false);
            std::cout << "Integrating\n";
            std::cout << "TImestep is " << dt << "\n";
            redraw();
            double min_l = 1e5;
            double max_l = 0.0;
            double l;
            for (Edge e : mesh->edges())
            {
                l = geometry->edgeLength(e);
                if (l < min_l)
                    min_l = l;
                if (l > max_l)
                    max_l = l;
            }
            std::cout << "The min edge length is " << min_l << " and the max edge length is " << max_l << "\n";
        }
    }
    ImGui::InputInt("Remeshing operations", &remeshing_ops);
}

int main(int argc, char **argv)
{

    bool Eigenvals = true;
    std::fstream JsonFile;

    JsonFile.open(argv[1], std::ios::in);
    int Nsim = std::stoi(argv[2]);
    json Data = json::parse(JsonFile);
    std::string first_dir;

    bool loaded = false;
    if (Data.contains("Subfolder"))
    {
        loaded = true;
        std::string subfolder = Data["Subfolder"];
        std::cout << "Loaded the subfolder " << subfolder << "\n";
        first_dir = Data["first_dir"];
        basic_name = first_dir + subfolder;
        std::cout << "THe basic name is " << basic_name << "\n";
        JsonFile.close();
        std::cout << "Closed the json\n";
        std::cout << "The candidate for file is " << basic_name + "Input_file.json" << "\n";
        JsonFile.open(basic_name + "/Input_file.json", std::ios::in);
        Data = json::parse(JsonFile);
    }

    Switch = "None";
    Switch_t = 0;

    bool finish_sim = false;
    if (Data.contains("finish_sim"))
    {
        finish_sim = Data["finish_sim"];
    }
    else
    {
        std::cout << "The simulation will finish when the timestep decreases\n";
    }

    std::string filepath = Data["init_file"];
    save_interval = Data["save_interval"];
    bool remesher = Data["remeshing"];

    if (Data.contains("Integration"))
    {
        Integration = Data["Integration"];
    }
    else
    {
        std::cout << "The integration method is not defined, using Gradient descent\n";
    }
    int Stored_info = 10;
    if (Data.contains("BFGS_saved_states"))
    {
        Stored_info = Data["BFGS_saved_states"];
    }

    // Switches handling
    std::vector<std::string> Switches(0);
    std::string Switch = "None";
    int Switch_t = 0;
    std::unordered_map<std::string, int> Switch_times_map;

    if (Data.contains("Switches"))
    {
        for (auto sw : Data["Switches"])
            Switches.push_back(sw);

        int switch_counter = 0;
        for (auto t : Data["Switch_times"])
        {
            Switch_times_map[Switches[switch_counter]] = t;
            switch_counter += 1;
        }
    }
    else
    {
        std::cout << "No switches in this run";
    }

    std::cout << "THe number of switches is " << Switches.size() << "\n";
    for (int i = 0; i < Switches.size(); i++)
    {
        std::cout << "The switch " << Switches[i] << " happens at time " << Switch_times_map[Switches[i]] << "\n";
    }

    if (Data.contains("remesh_every"))
        remesh_every = Data["remesh_every"];
    std::cout << "Remeshing every " << remesh_every << " steps " << std::endl;

    if (Data.contains("adapt_remesh"))
        adapt_remesh = Data["adapt_remesh"];

    Vector3 Recenter{0.0, 0.0, 0.0};
    if (Data.contains("Displacement"))
    {
        Recenter.x = Data["Displacement"][0];
        Recenter.y = Data["Displacement"][1];
        Recenter.z = Data["Displacement"][2];
        std::cout << "Displacing the membrane by " << Recenter << " \n";
    }

    double scale_factor = 1.0;
    if (Data.contains("rescale"))
    {
        scale_factor = Data["rescale"];
    }

    bool Saving_last_states = Data["saving_states"];
    size_t Final_t = Data["timesteps"];

    // Here i will load the geometry
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    // Recenter and rescale.

    if (Recenter.norm() > 0)
    {
        geometry->normalize(Recenter);
    }
    geometry->rescale(scale_factor);
    geometry->refreshQuantities();

    V_bar = geometry->totalVolume();
    A_bar = geometry->totalArea();

    // We will deal with the energies now

    for (auto Energy : Data["Energies"])
    {
        Energies.push_back(Energy["Name"]);
        Constants = Energy["constants"].get<std::vector<double>>();
        std::cout << "The constants for " << Energy["Name"] << " are ";
        for (size_t z = 0; z < Constants.size(); z++)
            std::cout << Constants[z] << " ";
        std::cout << " \n";

        // I want to add here something
        if (Energy["Name"] == "Volume_constraint")
        {
            Constraints.push_back("Volume");
            if (Constants[1] < 0)
            {
                Constants[1] = geometry->totalVolume();
                std::cout << "Setting the target volume to current volume \n";
            }
            V_bar = Constants[1];
        }
        if (Energy["Name"] == "Area_constraint")
        {
            Area_constraint = true;
            // std::cout << "\tThe area constraint is active\n";
            Constraints.push_back("Area");
            // The format for the area constraint is KA A_bar nu dA
            origA = geometry->totalArea();
            if (Constants[2] > 0)
            {
                double nu = Constants[2];
                A_bar = pow(36 * PI * V_bar * V_bar / (nu * nu), 1.0 / 3.0);
                Area = geometry->totalArea();
                std::cout << "The target Area is " << A_bar << " from nu " << nu << "\n";
                std::cout << "The current Area is " << Area << "\n";

                dA = Constants[3];
                if ((A_bar - Area) * dA < 0)
                    dA = -dA;
                minStepsReq = int(fabs((A_bar - origA) / dA));
                // First problem
                Constants[1] = Area + (dA / (fabs(dA))) * std::min(fabs(dA), fabs(A_bar - Area));
            }
            else
            {
                if (Constants[1] > 0)
                {
                    A_bar = Constants[1];
                    dA = 0;
                    minStepsReq = 0;
                }
                else
                {
                    dA = 0;
                    minStepsReq = 0;
                    A_bar = geometry->totalArea();
                    Constants[1] = A_bar;
                }
            }
            //

            std::cout << "The target area is " << A_bar << "\n";
            std::cout << "The current reduced volume is " << 3 * V_bar / (4 * PI * pow((Area / (4 * PI)), 1.5)) << "\n";
            // Ok but the thing is that we cannot do A_bar, because A_bar can be too different, so we need a A_bar_current
        }
        if (Energy["Name"] == "Membrane_tension")
        {
            // OK so
            double R0 = pow(0.75 * geometry->totalVolume() / PI, 1.0 / 3.0);
            // std::cout << "R0 is " << R0 << " \n";
            double A0 = 4 * PI * R0 * R0;
            Constants[1] = A0 * Constants[1];
        }
        Energy_constants.push_back(Constants);
        Constants.resize(0);
    }

    // We add the turning constraints
    Constraints.push_back("CMx");
    Constraints.push_back("CMy");
    Constraints.push_back("CMz");
    Constraints.push_back("Rx");
    Constraints.push_back("Ry");
    Constraints.push_back("Rz");

    Vector3 BPos;
    double radius;
    double interaction_str;
    std::string state;
    std::string interaction_mem;
    std::string Constraint;
    std::vector<double> Constraint_constants;
    Bead PBead;
    Interaction *PInteraction;

    int bead_counter = 0;
    for (auto Bead_data : Data["Beads"])
    {
        std::cout << "Adding a bead\n";
        std::cout << "THe bead is " << Bead_data["mem_inter"] << " \n";
        if (Bead_data.contains("gradient_order"))
        {
            Energies.push_back(Bead_data["gradient_order"]);
        }
        else
        {
            Energies.push_back("Bead");
        }

        if (Bead_data.contains("Constraint"))
        {
            Constraint = Bead_data["Constraint"];
            Constraint_constants = Bead_data["Constraint_constants"].get<std::vector<double>>();
        }
        else
        {
            Constraint = "None";
            Constraint_constants = {};
        }

        Energy_constants.push_back(Constants);

        BPos = Vector3({Bead_data["Pos"][0], Bead_data["Pos"][1], Bead_data["Pos"][2]});
        radius = Bead_data["radius"];
        state = Bead_data["state"];
        interaction_str = Bead_data["inter_str"];
        interaction_mem = Bead_data["mem_inter"];

        std::vector<double> Bead_params(0);
        Bead_params.push_back(interaction_str);
        Bead_params.push_back(radius);
        // Bead_params.push_back(2.0);

        if (Bead_data.contains("rc"))
        {
            Bead_params.push_back(Bead_data["rc"]);
        }
        else
        {

            if (Bead_data["mem_inter"] == "LJ")
            {
                Bead_params.push_back(radius * pow(2, 1.0 / 6.0));
            }
            else if (Bead_data["mem_inter"] == "Frenkel" || Bead_data["mem_inter"] == "Frenkel_Normal_nopush")
            {
                Bead_params.push_back(radius * 2.0);
            }
            else
            {
                Bead_params.push_back(-1);
            }
        }

        if (interaction_mem == "Frenkel")
        {

            if (Bead_data.contains("outside"))
            {
                Bead_params.push_back(Bead_data["outside"]);
            }
            else
            {
                Bead_params.push_back(1.0); // default outside
            }

            Interaction_container.push_back(std::move(make_unique<Frenkel>(mesh, geometry, Bead_params)));

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;
            std::cout << "Trivial assignments done\n";
            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            std::cout << "Assigned INter \n";
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            std::cout << "Assined bead of inter\n";
            Beads[bead_counter].Bead_id = bead_counter;

            std::cout << "The energy constants are ";
            for (size_t i = 0; i < Interaction_container[bead_counter].get()->Energy_constants.size(); i++)
            {
                std::cout << Interaction_container[bead_counter].get()->Energy_constants[i] << " ";
            }
            std::cout << "\n";
        }

        if (interaction_mem == "Frenkel_Normal_nopush")
        {

            if (Bead_data.contains("outside"))
            {
                Bead_params.push_back(Bead_data["outside"]);
            }
            else
            {
                Bead_params.push_back(1.0); // default outside
            }

            Interaction_container.push_back(std::move(make_unique<Frenkel_Normal>(mesh, geometry, Bead_params)));

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;
            std::cout << "Trivial assignments done\n";
            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            std::cout << "Assigned INter \n";
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            std::cout << "Assined bead of inter\n";
            Beads[bead_counter].Bead_id = bead_counter;

            std::cout << "The energy constants are ";
            for (size_t i = 0; i < Interaction_container[bead_counter].get()->Energy_constants.size(); i++)
            {
                std::cout << Interaction_container[bead_counter].get()->Energy_constants[i] << " ";
            }
            std::cout << "\n";
        }

        if (interaction_mem == "Linear")
        {

            if (Bead_data.contains("outside"))
            {
                Bead_params.push_back(Bead_data["outside"]);
            }
            else
            {
                Bead_params.push_back(1.0); // default outside
            }

            Interaction_container.push_back(std::move(make_unique<Linear>(mesh, geometry, Bead_params)));

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;
            std::cout << "Trivial assignments done\n";
            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            std::cout << "Assigned INter \n";
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            std::cout << "Assined bead of inter\n";
            Beads[bead_counter].Bead_id = bead_counter;

            std::cout << "The energy constants are ";
            for (size_t i = 0; i < Interaction_container[bead_counter].get()->Energy_constants.size(); i++)
            {
                std::cout << Interaction_container[bead_counter].get()->Energy_constants[i] << " ";
            }
            std::cout << "\n";
        }
        if (interaction_mem == "Linear_Normal")
        {

            if (Bead_data.contains("outside"))
            {
                Bead_params.push_back(Bead_data["outside"]);
            }
            else
            {
                Bead_params.push_back(1.0); // default outside
            }

            Interaction_container.push_back(std::move(make_unique<Linear_Normal>(mesh, geometry, Bead_params)));

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;
            std::cout << "Trivial assignments done\n";
            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            std::cout << "Assigned INter \n";
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            std::cout << "Assined bead of inter\n";
            Beads[bead_counter].Bead_id = bead_counter;

            std::cout << "The energy constants are ";
            for (size_t i = 0; i < Interaction_container[bead_counter].get()->Energy_constants.size(); i++)
            {
                std::cout << Interaction_container[bead_counter].get()->Energy_constants[i] << " ";
            }
            std::cout << "\n";
        }

        if (interaction_mem == "LJ")
        {

            if (Bead_data.contains("shift"))
            {
                Bead_params.push_back(Bead_data["shift"]);
            }
            else
            {
                Bead_params.push_back(0.0); // default shift
            }
            std::cout << "The bead params are " << Bead_params[0] << " " << Bead_params[1] << " " << Bead_params[2] << " " << Bead_params[3] << "\n";

            Interaction_container.push_back(std::move(make_unique<LJ>(mesh, geometry, Bead_params)));

            // std::cout<<"THe frenkel uptr points to"<< Interaction_container[bead_counter].get() <<"\n";

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;
            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            Beads[bead_counter].Bead_id = bead_counter;

            std::cout << "The energy at the cutfoff is" << Beads[bead_counter].Bead_I->E_r(Bead_params[2], Bead_params) << " and at more than the cutoff is " << Beads[bead_counter].Bead_I->E_r(1.5 * Bead_params[2], Bead_params) << "\n";

            std::cout << "The dE_r at the cutfoff is" << Beads[bead_counter].Bead_I->dE_r(Bead_params[2], Bead_params) << " and at more than the cutoff is " << Beads[bead_counter].Bead_I->dE_r(1.5 * Bead_params[2], Bead_params) << "\n";

            std::cout << "The ddE_r at the cutfoff is" << Beads[bead_counter].Bead_I->ddE_r(Bead_params[2], Bead_params) << " and at more than the cutoff is " << Beads[bead_counter].Bead_I->ddE_r(1.5 * Bead_params[2], Bead_params) << "\n";
        }
        if (interaction_mem == "One_over_r_x")
        {

            if (Bead_data.contains("outside"))
            {
                Bead_params.push_back(Bead_data["outside"]);
            }
            else
            {
                Bead_params.push_back(1.0); // default outside
            }
            Interaction_container.push_back(std::move(make_unique<One_over_r_Normal>(mesh, geometry, Bead_params)));

            std::cout << "THe  uptr points to" << Interaction_container[bead_counter].get() << "\n";

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;

            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            Beads[bead_counter].Bead_id = bead_counter;
        }

        if (interaction_mem == "One_over_r")
        {

            Interaction_container.push_back(std::move(make_unique<One_over_r>(mesh, geometry, Bead_params)));

            std::cout << "THe  uptr points to" << Interaction_container[bead_counter].get() << "\n";

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;

            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            Beads[bead_counter].Bead_id = bead_counter;
        }

        if (interaction_mem == "None")
        {
            std::cout << "Adding a no mem interaction\n";
            Interaction_container.push_back(std::move(make_unique<No_mem_Inter>()));

            Beads.push_back(Bead());
            Beads[bead_counter].mesh = mesh;
            Beads[bead_counter].geometry = geometry;
            Beads[bead_counter].Pos = BPos;
            Beads[bead_counter].strength = Bead_params[0];
            Beads[bead_counter].sigma = Bead_params[1];
            Beads[bead_counter].rc = Bead_params[2];
            Beads[bead_counter].interaction = interaction_mem;

            Beads[bead_counter].Bead_I = Interaction_container[bead_counter].get();
            Beads[bead_counter].Bead_I->Bead_1 = &Beads[bead_counter];
            Beads[bead_counter].Bead_id = bead_counter;
        }
        // BEA.interaction = interaction_mem;
        Beads[bead_counter].Bond_type = Bead_data["bonds"].get<vector<std::string>>();
        Beads[bead_counter].Interaction_constants_vector = Bead_data["bonds_constants"].get<std::vector<std::vector<double>>>();
        Beads[bead_counter].state = state;

        Beads[bead_counter].Constraint = Constraint;
        Beads[bead_counter].Constraint_constants = Constraint_constants;
        std::cout << "The bead has radius" << Beads[bead_counter].sigma << " cutoff of " << Beads[bead_counter].rc << " \n";

        // Lets add the manual movement here
        if (Beads[bead_counter].state == "manual")
        {
            // We need to set the velocity of the bead
            Beads[bead_counter].Velocity = Vector3({Bead_data["Velocity"][0], Bead_data["Velocity"][1], Bead_data["Velocity"][2]});
            if (Bead_data.contains("FinalPos"))
            {
                Beads[bead_counter].FinalPos = Vector3({Bead_data["FinalPos"][0], Bead_data["FinalPos"][1], Bead_data["FinalPos"][2]});
            }
        }

        bead_counter += 1;
    }
    std::cout << "The bead counter is " << bead_counter << "\n";

    for (int bi = 0; bi < Beads.size(); bi++)
    {
        std::cout << "The bead " << bi << " points to " << Beads[bi].Bead_I->Bead_1 << " and state " << Beads[bi].state << " \n";
    }

    // M3DG.Sim_handler->Calculate_ener

    std::vector<int> BeadBonds(0);
    size_t counter = 0;
    for (auto Bead_data : Data["Beads"])
    {
        // I am iterating again but i only care about the bonds
        BeadBonds = Bead_data["Beads"].get<std::vector<int>>();
        for (size_t i = 0; i < BeadBonds.size(); i++)
        {
            Beads[counter].Beads.push_back(&Beads[BeadBonds[i]]);
            counter += 1;
        }
    }
    std::cout << "There are " << counter / 2 << " bonds \n";

    if (counter > 0)
        std::cout << "The bead with radius " << Beads[0].sigma << " is connected to the bead with radius " << Beads[0].Beads[0]->sigma << " \n";

    for (size_t i = 0; i < Beads.size(); i++)
    {
        std::cout << "The bead has " << Beads[i].Bond_type.size() << " bonds and interaction " << Beads[i].interaction << " \n";
        std::cout << "The bead is at position " << Beads[i].Pos << "\n";
    }
    // Now the beads point to each other ()

    // Lets define our integrator and all its values

    int N_beads = Beads.size();
    for (int i = 0; i < N_beads; i++)
    {
        Beads[i].Total_beads = N_beads;
        Beads[i].Bead_id = i;
    }

    M3DG = Mem3DG(mesh, geometry);
    Sim_handler = E_Handler(mesh, geometry, Energies, Energy_constants);
    Og_Energy_constants = Energy_constants;
    Sim_handler.Trgt_vol = V_bar;
    Sim_handler.Trgt_area = A_bar;
    M3DG.recentering = Data["recentering"];
    M3DG.boundary = Data["boundary"];
    Sim_handler.boundary = Data["boundary"];

    for (size_t i = 0; i < Beads.size(); i++)
    {
        M3DG.Add_bead(&Beads[i]);
        Sim_handler.Add_Bead(&Beads[i]);
        // Beads[i].mesh = mesh;
        // Beads[i].geometry = geometry;
    }

    M3DG.Sim_handler = &Sim_handler;

    if (Data.contains("Field"))
    {
        M3DG.Field = Data["Field"];
        M3DG.Field_vals = Data["Field_vals"].get<std::vector<double>>();
        std::cout << "The field is " << M3DG.Field << " and the values are ";
        for (size_t i = 0; i < M3DG.Field_vals.size(); i++)
            std::cout << M3DG.Field_vals[i] << " ";
        std::cout << "\n";
    }
    else
    {
        M3DG.Field = "None";
    }

    if (Data.contains("backtrack"))
    {
        M3DG.backtrack = Data["backtrack"];
        if (Data.contains("Timestep"))
        {
            M3DG.timestep = Data["Timestep"];
        }
        else
        {
            M3DG.timestep = 1e-4;
        }
    }
    else
    {
        M3DG.backtrack = true;
    }

    if (Data.contains("momentum"))
    {
        M3DG.momentum = Data["momentum"];
        std::cout << "Momentum is " << M3DG.momentum << "\n";
    }

    if (Data.contains("remesher"))
    {
        Options.max_absolute_length = Data["remesher"]["size_max"];
        Options.min_absolute_length = Data["remesher"]["size_min"];
        Options.refine_angle = Data["remesher"]["refine_angle"];
        Options.aspect_min = Data["remesher"]["aspect_min"];
        Options.maxIterations = 1;
    }

    int saved_mesh_idx = 0;
    std::vector<Vector3> Bead_pos_saved(6);

    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();

    auto start_full = chrono::steady_clock::now();
    auto end_full = chrono::steady_clock::now();

    auto start_time_control = chrono::steady_clock::now();
    auto end_time_control = chrono::steady_clock::now();

    double remeshing_elapsed_time = 0;
    double integrate_elapsed_time = 0;
    double saving_mesh_time = 0;

    std::cout << "Current path is " << argv[0];

    std::cout << "The energy elements are \n";

    for (size_t z = 0; z < Energies.size(); z++)
    {
        std::cout << Energies[z] << " ";
    }
    std::cout << "\n";

    std::cout << "Minimum edge length allowed is " << Options.min_absolute_length << " muak\n";

    double avg_dih = 0;
    double max_dih = 0;
    double min_dih = 0.1;
    double dih;

    for (Edge e : mesh->edges())
    {
        dih = fabs(geometry->dihedralAngle(e.halfedge()));
        avg_dih += dih;
        if (dih > max_dih)
            max_dih = dih;
        if (dih < min_dih)
            min_dih = dih;
    }
    std::cout << "First checkThe average dihedral is" << avg_dih / mesh->nEdges() << " \n";
    std::cout << "The min dih is" << min_dih << " and the max dih is " << max_dih << " \n";
    avg_dih = 0.0;

    FaceData<double> F_sizings = M3DG.Face_sizings();
    VertexData<double> Sizings = M3DG.Vert_sizing(F_sizings);
    double max_sizing = 0.0;
    double min_sizing = 1e4;
    double sizing;

    for (Vertex v : mesh->vertices())
    {
        //
        sizing = Sizings[v];
        if (sizing > max_sizing)
            max_sizing = sizing;
        if (sizing < min_sizing)
            min_sizing = sizing;
    }
    std::cout << "My algorithm says\n";
    std::cout << "The max sizing is " << max_sizing << " and the min sizing is " << min_sizing << " \n";

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();

    // How do we create a name that makes sense
    // I can do something like E _ param_ param_ E_ param_param _Nsim_

    std::string Directory = "";
    std::stringstream stream;
    std::string s;
    bead_counter = 0;
    for (size_t z = 0; z < Energies.size(); z++)
    {
        Directory = Directory + Energies[z] + "_";

        if (Energies[z] == "Bead" || Energies[z] == "H1_Bead" || Energies[z] == "H2_Bead")
        {
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[bead_counter].sigma;
            Directory = Directory + "radius_" + stream.str() + "_";

            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[bead_counter].strength;
            Directory = Directory + "str_" + stream.str() + "_";

            if (Beads[bead_counter].Constraint == "Radial")
            {
                stream.str(std::string());
                stream << std::fixed << std::setprecision(4) << Beads[bead_counter].Constraint_constants[0];
                Directory = Directory + "theta_const_" + stream.str() + "_";
            }
            bead_counter += 1;
        }
        for (size_t j = 0; j < Energy_constants[z].size(); j++)
        {
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Energy_constants[z][j];

            Directory = Directory + stream.str() + "_";
            // stream << std::fixed << std::setprecision(2) << pi;
        }
    }
    if (Data.contains("Field"))
    {

        Directory = Directory + M3DG.Field + "_";
        for (size_t i = 0; i < M3DG.Field_vals.size(); i++)
        {
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << M3DG.Field_vals[i];
            Directory = Directory + stream.str() + "_";
        }
    }

    Count_remesh = false;
    if (Data.contains("Count_remesh"))
    {
        Count_remesh = Data["Count_remesh"];
        std::cout << "The remeshing counting is " << Count_remesh << "\n";
        Remeshing_count = std::ofstream(basic_name + "Remeshing_count.txt", std::ios_base::app);
        Remeshing_count << "#### timestep remeshing_operations nVertices nEdges nFaces \n";
        Remeshing_count.close();
    }

    bool bonds_exist = false;
    for (size_t z = 0; z < Beads.size(); z++)
    {
        if (Beads[z].Bond_type.size() > 0 && bonds_exist == false)
        {
            Directory = Directory + "Bonds_";
            bonds_exist = true;
        }
        for (size_t j = 0; j < Beads[z].Bond_type.size(); j++)
        {
            stream.str(std::string());
            stream << std::fixed << std::setprecision(4) << Beads[z].Interaction_constants_vector[j][0];
            Directory = Directory + Beads[z].Bond_type[j] + "_" + stream.str() + "_";
        }
    }

    if (Switch != "None")
    {
        Directory = Directory + "Switch_" + Switch + "_";
        stream.str(std::string());
        stream << std::fixed << std::setprecision(1) << Switch_t;
        Directory = Directory + "Switch_t_" + stream.str() + "_";
    }

    Directory = Directory + "Nsim_" + std::to_string(Nsim) + "/";

    std::cout << "Directory is " << Directory << " \n";

    //

    first_dir = Data["first_dir"];

    int status = mkdir(first_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    // std::string filename;
    if (loaded)
    {
        filename = basic_name + "Output_data.txt";
        Sim_data = std::ofstream(filename, std::ios_base::app);
        Sim_data.close();
    }
    else
    {

        basic_name = first_dir + Directory;

        status = mkdir(basic_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        std::cout << "\nIf this number is 0 the directory was created succesfully " << status << "\n";

        filename = basic_name + "Output_data.txt";

        Sim_data = std::ofstream(filename, std::ios_base::app);
        Sim_data << "time step Volume Area ";
        for (size_t i = 0; i < Energies.size(); i++)
        {
            Sim_data << Energies[i] << " ";
        }
        Sim_data << " Total_E grad_norm backtrackstep\n";
        Sim_data.close();
    }

    std::cout << "Here1\n";

    // std::vector<std::string> Bead_filenames;
    // std::ofstream Bead_datas;

    std::cout << "Here2\n";

    for (size_t i = 0; i < Beads.size(); i++)
    {
        Bead_filenames.push_back(basic_name + "Bead_" + std::to_string(i) + "_data.txt");
        Bead_datas = std::ofstream(Bead_filenames[i], std::ios_base::app);
        if (!Data.contains("Subfolder"))
        {
            Bead_datas << "####### This data is taken ever y" << save_interval << " steps just like the mesh radius is " << radius << " \n";
        }
        Bead_datas.close();
    }

    Count_remesh = false;
    std::ofstream Remeshing_count;
    if (Data.contains("Count_remesh"))
    {
        Count_remesh = Data["Count_remesh"];
        std::cout << "The remeshing counting is " << Count_remesh << "\n";
        Remeshing_count = std::ofstream(basic_name + "Remeshing_count.txt", std::ios_base::app);
        Remeshing_count << "#### timestep remeshing_operations nVertices nEdges nFaces \n";
        Remeshing_count.close();
    }

    std::cout << "Here3\n";

    Bead_filenames.push_back(basic_name + "Simulation_timings.txt");
    Bead_datas = std::ofstream(Bead_filenames[Beads.size()], std::ios_base::app);
    Bead_datas << "Remeshing_time Gradients Backtracking  Construction compute Solve \n";
    Bead_datas.close();

    std::cout << "Here4\n";

    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

    polyscope::TransparencyMode::Pretty;
    polyscope::init();
    TOTAL_ANGLE_DEFECT = geometry->totalAngleDefect();
    EULER_CHARACTERISTIC = geometry->eulerCharacteristic();

    std::cout << "Im here happy\n";

    DIR *FD;
    struct dirent *in_file;
    const char *target_dir = basic_name.c_str();

    FD = opendir(target_dir);
    if (FD == NULL)
    {
        fprintf(stderr, "Error : Failed to open target directory - %s\n", strerror(errno));
        return 1;
    }

    /* Reading object (files, directories ...) from the folder) */
    while ((in_file = readdir(FD)))
    {

        if (DT_REG == in_file->d_type)
        {
            if (in_file->d_name[0] == 'm')
            {
                // Ok SO THIS WILL GIVE ME
                std::vector<std::string> split_result = split(in_file->d_name, "_");
                // Ok now i want the second element
                int frame = std::stoi(split_result[1]);
                frames.push_back(frame);
                // mesh_names.push_back(in_file->d_name);
                // printf("%s file added\n", in_file->d_name);
            }
        }
    }

    std::sort(frames.begin(), frames.end());
    for (int i = 0; i < frames.size(); i++)
    {
        mesh_names.push_back("membrane_" + std::to_string(frames[i]) + ".obj");
        // std::cout << "The frame is " << frames[i] << "\n";
    }
    // return 1;

    /* Close the directory */
    closedir(FD);

    // OK i want to load the bead POSITIONS
    std::vector<std::pair<size_t, Vector3>> BeadPosPair;
    std::ifstream BeadReadData;
    size_t frame;
    Vector3 beadPos;
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
        // I want to do it for each one
        BeadPosPair.resize(0);
        // OK now i want to add the pairs
        BeadReadData = std::ifstream(basic_name + "Bead_" + std::to_string(bi) + "_data.txt");
        std::string line;
        while (std::getline(BeadReadData, line))
        {
            std::vector<std::string> splitted = split(line, " ");
            if (line[0] == '#')
                continue;
            // If not i need to use it
            frame = std::stoi(splitted[0]);
            beadPos = Vector3{std::stod(splitted[1]), std::stod(splitted[2]), std::stod(splitted[3])};
            BeadPosPair.push_back(std::make_pair(frame, beadPos));
        }
        std::sort(BeadPosPair.begin(), BeadPosPair.end(), deterministic_sort);
        // now we have a list of frames that is being sorter
        frameBeadPos.push_back(BeadPosPair);
    }

    // std::sort(mesh_names.begin(), mesh_names.end());

    // return 1;

    // polyscope::state::userCallback = Callback_curvatures;
    // polyscope::state::userCallback = Callback_anim;

    polyscope::state::userCallback = Callback_qts;

    geometry->requireVertexPositions();

    psMesh = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
    psMesh->setEdgeColor({0.0, 0.0, 0.0});
    psMesh->setTransparency(0.7);
    std::cout << "Registeredmesh\n";
    std::vector<glm::vec3> points;
    for (size_t i = 0; i < Beads.size(); i++)
    {
        points.push_back(glm::vec3(Beads[i].Pos.x, Beads[i].Pos.y, Beads[i].Pos.z));
        std::cout << "Points is " << points[i].x << " " << points[i].y << " " << points[i].z << "\n";
    }
    // std::cout<<"Points is "<< points[0][0]<<" "<< points[0][1] <<" "<< points[0][2] << "\n";
    if (Beads.size() > 0)
    {
        psCloud = polyscope::registerPointCloud("Great beads", points);
        psCloud->setPointRadius(Beads[0].sigma, false);
    }

    // std::cout<<"The radius would be " <<Beads[0].sigma*1.0 <<"\n";
    // psCloud->setPointRadius(Beads[0].sigma);
    // }

    std::cout << "Here7\n";

    psMesh->setSmoothShade(false);
    psMesh->setSurfaceColor({0.11, 0.388, 0.89});

    // i JUST WANT SOMETHING ELSE HERE
    double lengthScale = geometry->meanEdgeLength();
    vertexRadius = lengthScale * 0.02;
    edgeRadius = lengthScale * 0.01;

    std::cout << "Calling the mesh\n";
    polyscope::show();

    return EXIT_SUCCESS;
}
