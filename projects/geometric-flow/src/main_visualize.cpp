// GEOMETRIC FLOW

// #include <stdlib.h>
#include <unistd.h>

// #include <omp.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

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
#include "polyscope/view.h"
#include "glm/gtc/matrix_transform.hpp"

#include "args/args.hxx"
#include "imgui.h"

#include "Mem-3dg.h"
#include "Beads.h"
#include "Interaction.h"
#include "Energy_Handler.h"
#include "SimConfig.h"
#include "math.h"


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

// Everything loaded from the input file lives in sim (SimConfig.h). The
// viewer callbacks use the short names below, which alias into it.
Simulation sim;
ManifoldSurfaceMesh *&mesh = sim.mesh;
VertexPositionGeometry *&geometry = sim.geometry;
Mem3DG &M3DG = sim.M3DG;
E_Handler &Sim_handler = sim.Sim_handler;
std::vector<Bead> &Beads = sim.Beads;
std::vector<std::string> &Energies = sim.Energies;
std::vector<std::vector<double>> &Energy_constants = sim.Energy_constants;
RemeshOptions &Options = sim.Options;
double &V_bar = sim.V_bar;
double &A_bar = sim.A_bar;
double &dA = sim.dA;

// Used when loading saved frames
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;

polyscope::PointCloud *psCloud;
SimplePolygonMesh Saved_mesh;

// Polyscope has no plane primitive, so a bead with a Plane_adhesion
// interaction is drawn as a gray opaque quad of side PLANE_SIZE centred on the
// bead and perpendicular to the plane normal, and its sphere is hidden.
const double PLANE_SIZE = 10.0;

Plane_Adhesion *planeInteraction(size_t bi)
{
    return dynamic_cast<Plane_Adhesion *>(Beads[bi].Bead_I);
}

void updatePlanes()
{
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
        Plane_Adhesion *plane = planeInteraction(bi);
        if (plane == nullptr)
            continue;

        Vector3 n = plane->normal();
        Vector3 axis = fabs(n.x) < 0.9 ? Vector3{1.0, 0.0, 0.0} : Vector3{0.0, 1.0, 0.0};
        Vector3 u = unit(cross(n, axis)) * (0.5 * PLANE_SIZE);
        Vector3 v = unit(cross(n, u)) * (0.5 * PLANE_SIZE);
        const Vector3 &c = Beads[bi].Pos;

        // u, v, n right handed: the corners go counterclockwise seen from +n
        std::vector<Vector3> corners = {c - u - v, c + u - v, c + u + v, c - u + v};
        std::vector<std::vector<size_t>> quad = {{0, 1, 2, 3}};
        polyscope::SurfaceMesh *psPlane = polyscope::registerSurfaceMesh("Plane bead " + std::to_string(bi), corners, quad);
        psPlane->setSurfaceColor({0.5, 0.5, 0.5});
        psPlane->setTransparency(1.0);
        psPlane->setSmoothShade(false);
    }
}

// All beads go in the cloud (so bead indices stay the same); the ones drawn
// as a plane get radius 0
void registerBeadCloud(const std::vector<Vector3> &positions)
{
    psCloud = polyscope::registerPointCloud("Great beads", positions);
    psCloud->setPointRadius(Beads[0].sigma, false); // TODO adjust o different sizes

    std::vector<double> radii(Beads.size(), Beads[0].sigma);
    bool anyPlane = false;
    for (size_t bi = 0; bi < Beads.size(); bi++)
        if (planeInteraction(bi) != nullptr)
        {
            radii[bi] = 0.0;
            anyPlane = true;
        }
    if (anyPlane)
    {
        psCloud->addScalarQuantity("bead radius", radii);
        psCloud->setPointRadiusQuantity("bead radius", false);
    }
    updatePlanes();
}

// The plane is much bigger than the membrane, so the home view (which frames
// everything) would shrink the membrane. With a plane bead the camera looks at
// the membrane only, from the side of the plane the membrane is on, with the
// plane normal pointing up.
void frameCameraOnMembrane()
{
    Plane_Adhesion *plane = nullptr;
    for (size_t bi = 0; bi < Beads.size() && plane == nullptr; bi++)
        plane = planeInteraction(bi);
    if (plane == nullptr)
        return;

    geometry->requireVertexPositions();
    Vector3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
    for (Vertex v : mesh->vertices())
    {
        const Vector3 &p = geometry->inputVertexPositions[v];
        lo = Vector3{std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = Vector3{std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    Vector3 extent = hi - lo;
    Vector3 mid = 0.5 * (lo + hi);
    double size = std::max(extent.x, std::max(extent.y, extent.z));

    const Vector3 n = plane->normal();
    Vector3 side = fabs(n.y) < 0.9 ? Vector3{0.0, -1.0, 0.0} : Vector3{1.0, 0.0, 0.0};
    side = unit(side - dot(side, n) * n);
    const double elevation = 0.45; // radians above the plane
    Vector3 dir = cos(elevation) * side + sin(elevation) * n;

    Vector3 eye = mid + 2.2 * size * dir;
    polyscope::view::viewMat = glm::lookAt(glm::vec3(eye.x, eye.y, eye.z), glm::vec3(mid.x, mid.y, mid.z),
                                           glm::vec3(n.x, n.y, n.z));
}

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
int polish_left = 0; // gradient descent steps left in the current polish (PolishParams)
int n_polish = 0;
int remesh_pause_until = 0; // no remeshing before this integration step (after a rollback)
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


std::string basic_name;
int Save_slot = 0;
int Load_slot = 0;

std::vector<std::vector<double>> Og_Energy_constants(0);
std::vector<double> Constants(0);
std::ofstream Sim_data;
std::vector<std::string> Bead_filenames;
std::ofstream Bead_datas;
// std::vector<Be>
std::vector<std::string> Constraints(0);
std::vector<std::string> bonds;
std::vector<std::vector<double>> constants;
// std::vector<Interaction *> Bead_Interactions(0);

std::vector<std::string> mesh_names(0);
std::vector<int> frames(0);
std::vector<std::vector<std::pair<size_t, Vector3>>> frameBeadPos;
std::vector<Vector3> CurrentBeadPos;
std::vector<Vector3> CurrentBeadForce;

polyscope::SurfaceMesh *psMesh;

// Light blue of the paper renders. showSelected paints the unselected faces with
// it too, so a click on the membrane keeps the colour.
const std::array<double, 3> MEMBRANE_COLOR = {0.25, 0.55, 0.87};

// Polyscope forgets the look of a structure when it is registered again, and the
// membrane is registered after every step, touch and frame change, so the look
// is set here and not once at startup.
polyscope::SurfaceMesh *registerMembrane()
{
    polyscope::SurfaceMesh *ps = polyscope::registerSurfaceMesh("MyMesh", geometry->vertexPositions, mesh->getFaceVertexList());
    ps->setSmoothShade(false);
    ps->setSurfaceColor({MEMBRANE_COLOR[0], MEMBRANE_COLOR[1], MEMBRANE_COLOR[2]});
    ps->setEdgeColor({0.0, 0.0, 0.0});
    ps->setEdgeWidth(1.0);
    return ps;
}
double TOTAL_ANGLE_DEFECT;
size_t EULER_CHARACTERISTIC;

std::map<int, std::vector<std::array<double, 3>>> sColors;
std::map<int, std::vector<std::array<double, 3>>> sColors2;

polyscope::SurfaceVertexColorQuantity *vertexColors;
polyscope::SurfaceFaceColorQuantity *faceColors;
double vertexRadius;
double edgeRadius;

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
    std::vector<Vector3> edgePos;
    std::vector<std::array<size_t, 2>> edgeInd;
    for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
         it != polyscope::state::subset.edges.end(); ++it)
    {
        Edge e = mesh->edge(*it);
        // std::cout << "Showing edge " << *it << "\n";
        // std::cout << "Connected by " << e.firstVertex() << " and " << e.secondVertex() << "\n";
        // std::cout<<"The vertices are" <<
        edgePos.push_back(geometry->vertexPositions[e.firstVertex()]);
        edgePos.push_back(geometry->vertexPositions[e.secondVertex()]);
        size_t i = edgeInd.size();
        edgeInd.push_back({2 * i, 2 * i + 1});
    }
    polyscope::SurfaceGraphQuantity *showEdges = psMesh->addSurfaceGraphQuantity("selected edges", edgePos, edgeInd);
    showEdges->setEnabled(true);
    showEdges->setRadius(edgeRadius);
    showEdges->setColor(ORANGE_VEC);

    // Show selected faces.
    std::vector<std::array<double, 3>> faceColors(mesh->nFaces());
    for (size_t i = 0; i < mesh->nFaces(); i++)
    {
        faceColors[i] = MEMBRANE_COLOR;
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
    {
        psCloud->updatePointPositions(BeadPositions);
        updatePlanes();
    }
    polyscope::requestRedraw();
}

void Save_mesh(std::string basic_name, size_t current_t)
{
    Save_mesh(mesh, geometry, basic_name, current_t);
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
    if (Energy == "Coverage")
        return Sim_handler.F_Coverage(energyConstants);
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
        // After a remesh or a history reset too, as main_cluster does: a
        // remesh adds vertices that would otherwise have no normal
        if (First_iter || M3DG.BFGS_iter == 0)
        {
            std::cout << " Updating vertex normals for BFGS-Normal\n";
            M3DG.Sim_handler->update_vertex_normals();
        }
        if (First_iter)
        {
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
            Jt = Jt.block(0, 0, Sim_handler.mesh->nVertices(), 1).eval(); // eval: the block aliases Jt
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
    Jt = Jt.block(0, 0, Sim_handler.mesh->nVertices(), 2).eval(); // eval: the block aliases Jt
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

int currFrame = 0;
bool autoPlaying = false;

// The viewer shows one run folder (basic_name): its Input_file.json built the
// simulation in sim, and its membrane_<t>.obj and Bead_<i>_data.txt are the
// frames the slider goes through.

// Reads the frames saved in basic_name into frames, mesh_names and frameBeadPos
// (one row per saved frame for each bead). False, with a message, if the
// folder cannot be read.
bool read_frames()
{
    DIR *FD = opendir(basic_name.c_str());
    if (FD == NULL)
    {
        fprintf(stderr, "Error : Failed to open target directory - %s\n", strerror(errno));
        return false;
    }
    frames.resize(0);
    struct dirent *in_file;
    while ((in_file = readdir(FD)))
    {
        if (DT_REG == in_file->d_type && in_file->d_name[0] == 'm')
        {
            // membrane_<t>.obj: the frame is the second element
            std::vector<std::string> split_result = split(in_file->d_name, "_");
            if (split_result.size() > 1)
                frames.push_back(std::stoi(split_result[1]));
        }
    }
    closedir(FD);

    std::sort(frames.begin(), frames.end());
    mesh_names.resize(0);
    for (size_t i = 0; i < frames.size(); i++)
        mesh_names.push_back("membrane_" + std::to_string(frames[i]) + ".obj");

    frameBeadPos.resize(0);
    std::vector<std::pair<size_t, Vector3>> BeadPosPair;
    for (size_t bi = 0; bi < Beads.size(); bi++)
    {
        BeadPosPair.resize(0);
        std::ifstream BeadReadData(basic_name + "Bead_" + std::to_string(bi) + "_data.txt");
        std::string line;
        while (std::getline(BeadReadData, line))
        {
            if (line[0] == '#')
                continue;
            std::vector<std::string> splitted = split(line, " ");
            size_t frame = std::stoi(splitted[0]);
            Vector3 beadPos = Vector3{std::stod(splitted[1]), std::stod(splitted[2]), std::stod(splitted[3])};
            BeadPosPair.push_back(std::make_pair(frame, beadPos));
        }
        std::sort(BeadPosPair.begin(), BeadPosPair.end(), deterministic_sort);
        frameBeadPos.push_back(BeadPosPair);
    }
    return true;
}

// Sets the viewer up for the run held by sim: the settings it shows, the output
// files and the frames. The run folder is the one of the config
// (loaded_from_subfolder) or a new numbered one in first_dir. Throws
// std::runtime_error if that folder cannot be made, or, with check_frames, if it
// has no frame to show (nothing is written to it then).
void prepare_run(bool check_frames = false)
{
    const SimConfig &cfg = sim.cfg;
    const bool loaded = cfg.loaded_from_subfolder;

    Switch = "None";
    Switch_t = 0;
    Integration = cfg.integration;
    save_interval = cfg.save_interval;
    remesh_every = cfg.remesh_every;
    adapt_remesh = cfg.adapt_remesh;
    Count_remesh = cfg.count_remesh;
    Og_Energy_constants = Energy_constants;

    integrating = false;
    integration_counter = 0;
    continueLastCall = false;
    savedSteps = 0;
    last_remesh = 0;
    polish_left = 0;
    n_polish = 0;
    remesh_pause_until = 0;
    first_newton = true;

    // Constraints and area schedule used by the Newton tools and updateTargetArea
    Constraints.clear();
    Area_constraint = false;
    minStepsReq = 0;
    origA = geometry->totalArea();
    for (size_t i = 0; i < Energies.size(); i++)
    {
        if (Energies[i] == "Volume_constraint")
            Constraints.push_back("Volume");
        if (Energies[i] == "Area_constraint")
        {
            Area_constraint = true;
            Constraints.push_back("Area");
            minStepsReq = fabs(dA) > 0 ? size_t(fabs((A_bar - origA) / dA)) : 0;
        }
    }
    for (const char *name : {"CMx", "CMy", "CMz", "Rx", "Ry", "Rz"})
        Constraints.push_back(name);

    ORIG_VPOS = geometry->inputVertexPositions;
    CoM = geometry->centerOfMass();

    if (loaded)
    {
        basic_name = cfg.subfolder_dir;
    }
    else
    {
        basic_name = make_numbered_dir(cfg.first_dir);
        record_run(cfg.first_dir, basic_name, build_descriptive_name(sim), cfg.source_path);
    }
    std::cout << "The output directory is " << basic_name << "\n";
    CurrentBeadPos.resize(0);
    CurrentBeadForce.resize(0);
    if (!read_frames())
        throw std::runtime_error("Could not read the frames in " + basic_name);
    if (check_frames)
    {
        if (frames.empty())
            throw std::runtime_error("there is no membrane_<t>.obj in " + basic_name);
        for (size_t bi = 0; bi < Beads.size(); bi++)
            if (frameBeadPos[bi].size() < frames.size())
                throw std::runtime_error("Bead_" + std::to_string(bi) + "_data.txt has " + std::to_string(frameBeadPos[bi].size()) +
                                         " rows for " + std::to_string(frames.size()) + " frames in " + basic_name);
    }

    filename = basic_name + "Output_data.txt";
    Bead_filenames = open_output_files(sim, basic_name, loaded);

    // The header goes in once, not at every opening of the run
    const std::string timings_file = basic_name + "Simulation_timings.txt";
    Bead_filenames.push_back(timings_file);
    bool timings_empty = true;
    {
        std::ifstream probe(timings_file);
        timings_empty = !probe.is_open() || probe.peek() == std::ifstream::traits_type::eof();
    }
    Bead_datas = std::ofstream(timings_file, std::ios_base::app);
    if (timings_empty)
        Bead_datas << "Remeshing_time Gradients Backtracking  Construction compute Solve \n";
    Bead_datas.close();

    TOTAL_ANGLE_DEFECT = geometry->totalAngleDefect();
    EULER_CHARACTERISTIC = geometry->eulerCharacteristic();
}

// The run folders next to basic_name, to switch between similar simulations
// ("12/" -> "../Results/Name/12/" and its siblings).
struct RunFolder
{
    std::string name;        // without the trailing '/'
    std::string description; // from runs_index.txt, empty if not recorded
};
std::vector<RunFolder> run_folders;
std::string run_switch_error;

// With --load_chain the folder of a chain (with trailing '/'): its link*/<n> run folders are the runs
// of the menu, named "link03_s4.30/1". Empty: the runs are the siblings of the current one.
std::string chain_root;

// "a/b/12/" -> "a/b/"
std::string run_parent(std::string run_dir)
{
    while (run_dir.size() > 1 && run_dir.back() == '/')
        run_dir.pop_back();
    size_t slash = run_dir.find_last_of('/');
    return slash == std::string::npos ? "" : run_dir.substr(0, slash + 1);
}

// Where the names of the menu start: the chain folder, or the parent of the run
std::string run_menu_root(const std::string &run_dir)
{
    return chain_root.empty() ? run_parent(run_dir) : chain_root;
}

// "a/b/12/" -> "12" (in a chain "a/chain/link03_s4.30/1/" -> "link03_s4.30/1")
std::string run_name(const std::string &run_dir)
{
    std::string name = run_dir.substr(run_menu_root(run_dir).size());
    while (!name.empty() && name.back() == '/')
        name.pop_back();
    return name;
}

// Numbered folders first and in numeric order, then the others alphabetically
bool run_name_less(const std::string &a, const std::string &b)
{
    auto numeric = [](const std::string &s)
    { return !s.empty() && s.find_first_not_of("0123456789") == std::string::npos; };
    auto no_zeros = [](const std::string &s)
    {
        size_t k = s.find_first_not_of('0');
        return k == std::string::npos ? std::string("0") : s.substr(k);
    };
    const bool na = numeric(a), nb = numeric(b);
    if (na != nb)
        return na;
    if (na)
    {
        const std::string x = no_zeros(a), y = no_zeros(b);
        if (x.size() != y.size())
            return x.size() < y.size();
        if (x != y)
            return x < y;
    }
    return a < b;
}

// The runs of a chain in the order they ran: by link number, a waypoint (link02_s5.00_w1) before the
// link it leads to, then the runs inside the link
bool chain_name_less(const std::string &a, const std::string &b)
{
    auto number_after = [](const std::string &name, const std::string &tag)
    {
        const size_t k = name.find(tag);
        return k == std::string::npos ? std::numeric_limits<int>::max() : atoi(name.c_str() + k + tag.size());
    };
    const std::string link_a = a.substr(0, a.find('/')), link_b = b.substr(0, b.find('/'));
    for (const char *tag : {"link", "_w"})
    {
        const int x = number_after(link_a, tag), y = number_after(link_b, tag);
        if (x != y)
            return x < y;
    }
    if (link_a != link_b)
        return link_a < link_b;
    return run_name_less(a.substr(std::min(link_a.size() + 1, a.size())), b.substr(std::min(link_b.size() + 1, b.size())));
}

// Appends the sub-directories of parent (with trailing '/', or empty for the current
// directory) that hold an Input_file.json, as prefix + name
void add_run_folders(const std::string &parent, const std::string &prefix)
{
    std::map<std::string, std::string> descriptions;
    std::ifstream index(parent + "runs_index.txt");
    std::string line;
    while (std::getline(index, line))
    {
        std::vector<std::string> columns = split(line, "\t");
        if (columns.size() < 2)
            continue;
        std::string folder = columns[0];
        while (!folder.empty() && folder.back() == '/')
            folder.pop_back();
        descriptions[folder] = columns[1];
    }

    DIR *dir = opendir(parent.empty() ? "." : parent.c_str());
    if (dir == NULL)
        return;
    while (struct dirent *entry = readdir(dir))
    {
        const std::string name = entry->d_name;
        if (name[0] == '.')
            continue;
        struct stat info;
        if (stat((parent + name + "/Input_file.json").c_str(), &info) != 0)
            continue;
        run_folders.push_back({prefix + name, descriptions[name]});
    }
    closedir(dir);
}

// The runs of parent (with trailing '/', or empty for the current directory): its sub-directories that
// hold an Input_file.json or, for a chain, the ones inside its link folders ("link03_s4.30/1")
void scan_run_folders_in(const std::string &parent, bool chain = false)
{
    run_folders.clear();
    if (!chain)
        add_run_folders(parent, "");
    else if (DIR *dir = opendir(parent.c_str()))
    {
        std::vector<std::string> links;
        while (struct dirent *entry = readdir(dir))
        {
            const std::string name = entry->d_name;
            if (name[0] != '.')
                links.push_back(name);
        }
        closedir(dir);
        for (const std::string &link : links)
            add_run_folders(parent + link + "/", link + "/");
    }
    std::sort(run_folders.begin(), run_folders.end(),
              [chain](const RunFolder &a, const RunFolder &b)
              { return chain ? chain_name_less(a.name, b.name) : run_name_less(a.name, b.name); });
}

// The folders next to the run being shown (the runs of the chain with --load_chain)
void scan_run_folders()
{
    scan_run_folders_in(run_menu_root(basic_name), !chain_root.empty());
}

// The config of the run in run_dir (with trailing '/'), as load_config gives it
// for a file with "Subfolder". Throws if it cannot be read.
SimConfig run_config(const std::string &run_dir)
{
    SimConfig cfg = load_config(run_dir + "Input_file.json");
    cfg.loaded_from_subfolder = true;
    cfg.subfolder_dir = run_dir;

    // A run folder that was moved or copied (the viewer_runs of Analyze_*_chains.py, the links of a chain
    // from the cluster) still names its mesh by the old absolute path: the first frame of the folder, the
    // state the run started from, stands in for it
    struct stat info;
    if (stat(cfg.init_file.c_str(), &info) != 0 || !S_ISREG(info.st_mode))
    {
        const std::string local = run_dir + "membrane_0.obj";
        if (stat(local.c_str(), &info) == 0 && S_ISREG(info.st_mode))
        {
            std::cout << "init_file " << cfg.init_file << " is missing, using " << local << "\n";
            cfg.init_file = local;
        }
    }
    return cfg;
}

// Drops the simulation and everything of the viewer that depends on it, then
// builds the one of cfg. Throws if that fails, leaving the viewer without a run.
void open_run(const SimConfig &cfg)
{
    // Whatever holds MeshData of the mesh goes before the mesh
    ORIG_VPOS = VertexData<Vector3>();
    Saved_vertex_positions = VertexData<Vector3>();
    Vert_sizings = VertexData<double>();
    Face_sizings = FaceData<double>();
    Edge_sizings = EdgeData<double>();
    Gradient_vertex = VertexData<Vector3>();

    // The drawings of the beads of the old run (the mesh is registered again by the frame load)
    for (size_t bi = 0; bi < Beads.size(); bi++)
        polyscope::removeStructure("Plane bead " + std::to_string(bi), false);
    polyscope::removeStructure("Great beads", false);
    psCloud = nullptr;

    sim.reset();
    build_simulation(cfg, sim);
    prepare_run(true);
}

// Show the run folder `name` next to the current one, at the same frame (or the
// last one it has). If it cannot be opened the run we came from is kept (opened
// again if it was already dropped) and run_switch_error says why. The caller
// reloads the frame afterwards.
void switch_run(const std::string &name)
{
    const std::string previous = basic_name;
    const std::string target = run_menu_root(previous) + name + "/";
    run_switch_error.clear();

    SimConfig cfg;
    try
    {
        cfg = run_config(target);
    }
    catch (const std::exception &e)
    {
        run_switch_error = "Could not open " + target + ": " + e.what();
        std::cerr << run_switch_error << "\n";
        return;
    }
    try
    {
        open_run(cfg);
    }
    catch (const std::exception &e)
    {
        run_switch_error = "Could not open " + target + ": " + e.what();
        std::cerr << run_switch_error << "\n";
        try
        {
            open_run(run_config(previous));
        }
        catch (const std::exception &e2)
        {
            std::cerr << "Could not go back to " << previous << ": " << e2.what() << "\n";
            std::exit(EXIT_FAILURE);
        }
    }
    currFrame = std::min(currFrame, int(frames.size()) - 1);
    scan_run_folders();
}

static const ImVec4 pressColor = ImColor::HSV(1. / 7.0f, 0.6f, 0.6f); // gold
static const ImVec4 releaseColor{0.35f, 0.61f, 0.49f, 0.62f};         // default green

static ImVec4 off_currColor = pressColor;
static ImVec4 V_sizing_currColor = releaseColor;
static ImVec4 F_sizing_currColor = releaseColor;

std::vector<ImVec4 *> sState = {&off_currColor, &V_sizing_currColor, &F_sizing_currColor};

// The run folder selector at the top of the window: a list of the folders next
// to the current one, previous/next and a refresh. True if the run was switched,
// so the caller loads the frame.
bool run_folder_menu()
{
    const std::string current = run_name(basic_name);
    int current_index = -1;
    for (size_t i = 0; i < run_folders.size(); i++)
        if (run_folders[i].name == current)
            current_index = int(i);

    // The label first and a narrow combo, so the row fits in the panel
    std::string target;
    ImGui::AlignTextToFramePadding();
    ImGui::Text("Run folder");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(ImGui::GetFontSize() * (chain_root.empty() ? 8.0f : 13.0f)); // "link02_s5.00_w1/1"
    if (ImGui::BeginCombo("##run_folder", current.c_str()))
    {
        for (size_t i = 0; i < run_folders.size(); i++)
        {
            const bool selected = int(i) == current_index;
            if (ImGui::Selectable(run_folders[i].name.c_str(), selected) && !selected)
                target = run_folders[i].name;
            if (selected)
                ImGui::SetItemDefaultFocus();
            if (ImGui::IsItemHovered() && !run_folders[i].description.empty())
            {
                ImGui::BeginTooltip();
                ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
                ImGui::TextUnformatted(run_folders[i].description.c_str());
                ImGui::PopTextWrapPos();
                ImGui::EndTooltip();
            }
        }
        ImGui::EndCombo();
    }
    ImGui::SameLine();
    if (ImGui::Button("<##run") && current_index > 0)
        target = run_folders[current_index - 1].name;
    ImGui::SameLine();
    if (ImGui::Button(">##run") && current_index >= 0 && current_index + 1 < int(run_folders.size()))
        target = run_folders[current_index + 1].name;
    ImGui::SameLine();
    if (ImGui::Button("Refresh##run"))
        scan_run_folders();
    if (current_index >= 0)
        ImGui::Text("Run %d of %zu", current_index + 1, run_folders.size());
    if (!run_switch_error.empty())
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "%s", run_switch_error.c_str());

    if (target.empty())
        return false;
    switch_run(target);
    return true;
}

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
        psMesh = registerMembrane();
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

// ---- Look of the right panel ------------------------------------------------
// One blue accent (the blue of the membrane) for the main menus, buttons, sliders
// and checkboxes, all with a thin outline. Set once after polyscope::init (which
// sets its own style).
void applyMenuStyle()
{
    ImGuiStyle &st = ImGui::GetStyle();
    st.FrameRounding = 4.0f;
    st.GrabRounding = 4.0f;
    st.PopupRounding = 4.0f;
    st.ScrollbarRounding = 6.0f;
    st.FrameBorderSize = 1.0f; // the outline of buttons, headers, sliders and inputs
    st.FramePadding = ImVec2(6.0f, 3.5f);
    st.ItemSpacing = ImVec2(8.0f, 6.0f);
    st.IndentSpacing = 14.0f;
    st.GrabMinSize = 12.0f;

    const ImVec4 accent(0.20f, 0.45f, 0.80f, 1.0f);
    const ImVec4 accentHover(0.27f, 0.55f, 0.92f, 1.0f);
    const ImVec4 accentActive(0.33f, 0.62f, 1.00f, 1.0f);
    const ImVec4 outline(0.40f, 0.65f, 0.98f, 0.85f);
    ImVec4 *c = st.Colors;
    c[ImGuiCol_Border] = outline;
    c[ImGuiCol_Button] = ImVec4(accent.x, accent.y, accent.z, 0.75f);
    c[ImGuiCol_ButtonHovered] = accentHover;
    c[ImGuiCol_ButtonActive] = accentActive;
    c[ImGuiCol_Header] = ImVec4(accent.x, accent.y, accent.z, 0.75f);
    c[ImGuiCol_HeaderHovered] = accentHover;
    c[ImGuiCol_HeaderActive] = accentActive;
    c[ImGuiCol_FrameBg] = ImVec4(0.12f, 0.20f, 0.32f, 0.70f);
    c[ImGuiCol_FrameBgHovered] = ImVec4(0.18f, 0.30f, 0.48f, 0.80f);
    c[ImGuiCol_FrameBgActive] = ImVec4(0.22f, 0.37f, 0.58f, 0.90f);
    c[ImGuiCol_CheckMark] = accentActive;
    c[ImGuiCol_SliderGrab] = accentHover;
    c[ImGuiCol_SliderGrabActive] = accentActive;
}

// ---- Background -------------------------------------------------------------
// Dark by default; the checkbox at the top of the right panel switches it.
bool darkBackground = true;

void applyBackground()
{
    if (darkBackground)
        polyscope::view::bgColor = {0.08f, 0.08f, 0.10f, 0.0f};
    else
        polyscope::view::bgColor = {1.0f, 1.0f, 1.0f, 0.0f};
    polyscope::requestRedraw();
}

// A main menu: a full-width blue header bar with an outline. Like TreeNodeEx:
// call ImGui::TreePop() if it returns true.
bool menuSection(const char *label)
{
    ImGui::Spacing();
    return ImGui::TreeNodeEx(label, ImGuiTreeNodeFlags_Framed | ImGuiTreeNodeFlags_SpanAvailWidth |
                                        ImGuiTreeNodeFlags_FramePadding);
}

// A group inside a main menu: a quiet grey bar. Same contract as menuSection.
bool subSection(const std::string &label)
{
    ImGui::PushStyleColor(ImGuiCol_Header, ImVec4(0.30f, 0.33f, 0.38f, 0.55f));
    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, ImVec4(0.38f, 0.42f, 0.48f, 0.70f));
    ImGui::PushStyleColor(ImGuiCol_HeaderActive, ImVec4(0.45f, 0.50f, 0.57f, 0.80f));
    const bool open = ImGui::TreeNodeEx(label.c_str(), ImGuiTreeNodeFlags_Framed | ImGuiTreeNodeFlags_SpanAvailWidth);
    ImGui::PopStyleColor(3);
    return open;
}

void Callback_qts()
{
    frame_counter++;
    if (frame_counter >= 60)
        frame_counter = frame_counter % 60;
    ImGui::Text("The current frame is %i", frame_counter);
    if (ImGui::Checkbox("Dark background", &darkBackground))
        applyBackground();

    bool updateFrameData = false;
    bool currframeUpdated = run_folder_menu();
    bool integratorUpdate = false;
    bool updateFrames = false;
    if (ImGui::Button("Update frames"))
    {
        updateFrames = true;
    }

    if ((integration_counter == integration_steps && integrating) || updateFrames)
    {
        integrating = false;
        integratorUpdate = true;

        std::cout << "Done integrating, we will update the files available.\n";
        if (!read_frames())
            return;
        std::cout << "Filenames updated\n";
        if (!updateFrames)
            currFrame = currFrame + savedSteps;
        updateFrames = false;
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

        psMesh = registerMembrane();
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
            registerBeadCloud(CurrentBeadPos);
        }
    }

    ImGuiTreeNodeFlags flag2 = ImGuiTreeNodeFlags_None;

    if (menuSection("Integration"))
    {
        // Ok so i am here and I already decided there are two things i need to do
        // One is the integration settings and one is the actual integration
        ImGui::Text(("The current integration method is " + Integration).c_str());
        if (ImGui::Button("Integrate!"))
        {
            integrating = true;
            savedSteps = 0;
            integration_counter = 0;
            remesh_pause_until = 0;
            // last_remesh = 0;
            current_t = frames[currFrame];
            std::cout << "We are starting from t = " << current_t << " \n";
        }
        ImGui::SameLine(150);
        ImGui::Checkbox("Continue last call", &continueLastCall);

        if (subSection("Options"))
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
        ImGui::TreePop();
    }
    // How often do we check the small angle? every 10?

    if (integrating && frame_counter % 2 == 0)
    {

        // Ok now i can integrate every second frame(hopefully this helps us)
        integration_counter += 1;
        bool Save = integration_counter % save_interval == 0;
        // I should save before running
        // std::cout << "We are at t = " << current_t + integration_counter << " \n";

        // Undo of this step's remeshing (RemeshRollbackParams), armed before
        // the mesh is first changed
        const RemeshRollbackParams &rollback = sim.cfg.remesh_rollback;
        const bool remesh_paused = integration_counter < remesh_pause_until;
        MeshSnapshot undo;
        bool undo_armed = false; // stays set until the step after the remesh is done
        double E_undo = 0.0;
        std::vector<double> terms_undo;
        auto arm_undo = [&]()
        {
            if (undo_armed || rollback.max_rise == 0.0)
                return;
            undo = take_snapshot(sim);
            terms_undo = energy_terms(sim, E_undo);
            undo_armed = true;
        };

        // i should add the small angle flag
        bool flagSmallAngle = false;
        if (integration_counter % 5 == 0 && !remesh_paused)
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
            const PolishParams &polish = sim.cfg.polish;
            if (flagSmallAngle && Integration == "BFGS-Normal" && polish.gd_steps > 0 &&
                (polish.max_cycles < 0 || n_polish < polish.max_cycles))
            {
                // Remesh (below) and relax with gradient descent before going back
                n_polish++;
                polish_left = polish.gd_steps;
                Integration = "Gradient_descent";
                std::cout << "Polish " << n_polish << ": slivers at t = " << current_t + integration_counter << "\n";
            }
            if (flagSmallAngle)
            {
                arm_undo();
                remeshSmallAngles(*mesh, *geometry, Options);
                M3DG.BFGS_iter = 0;
                mesh->compress();
                geometry->refreshQuantities(); // i guess yes
            }
            geometry->unrequireCornerAngles();
            // we do the small angle check
        }
        if (!remesh_paused &&
            (((integration_counter - last_remesh) >= remesh_every && remesh_every > 0) || flagSmallAngle || M3DG.remesh_flag))
        {
            M3DG.remesh_flag = false;
            last_remesh = integration_counter;
            arm_undo();
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
        if (undo_armed)
        {
            double E = 0.0;
            const std::vector<double> terms = energy_terms(sim, E);
            const int i = exploded_term(terms_undo, terms, rollback.max_rise);
            if (i >= 0)
            {
                // BFGS-Normal recomputes the normals since BFGS_iter is 0
                restore_snapshot(sim, undo);
                undo_armed = false;
                M3DG.BFGS_iter = 0;
                remesh_pause_until = integration_counter + rollback.wait + 1;
                std::cout << "Remesh at t = " << current_t + integration_counter << " raised " << Energies[i] << " from "
                          << terms_undo[i] << " to " << terms[i] << ", undone; no remeshing for " << rollback.wait
                          << " steps\n";
            }
        }
        if (Area_constraint && updateArea)
            updateTargetArea(current_t + integration_counter, origA, A_bar, dA, minStepsReq);
        double timestepIter = Integration_step(integration_counter + current_t, Save, integration_counter == 1 && !continueLastCall);
        if (polish_left > 0 && Integration != "Gradient_descent")
            polish_left = 0; // the integrator was changed in the menu
        if (polish_left > 0 && timestepIter > 0 && --polish_left == 0)
        {
            Integration = sim.cfg.polish.then;
            M3DG.BFGS_iter = 0;
            std::cout << "Polish " << n_polish << ": relaxed at t = " << current_t + integration_counter << ", back to "
                      << Integration << "\n";
        }

        if (Save)
        {
            savedSteps += 1;
            Save_mesh(basic_name, current_t + integration_counter);
        }
        if (timestepIter <= 0 && undo_armed)
        {
            // The step after a remesh failed: undo the remesh instead of stopping
            restore_snapshot(sim, undo);
            M3DG.BFGS_iter = 0;
            M3DG.remesh_flag = false;
            M3DG.step_failed = false;
            remesh_pause_until = integration_counter + rollback.wait + 1;
            std::cout << "The step after the remesh at t = " << current_t + integration_counter
                      << " failed, remesh undone; no remeshing for " << rollback.wait << " steps\n";
        }
        else if (timestepIter <= 0)
        {
            std::cout << "Sim is not progresssing we will stop integrating\n";
            integration_counter = integration_steps;
        }
    }

    if (menuSection("Remeshing"))
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
            psMesh = registerMembrane();
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
            psMesh = registerMembrane();
            polyscope::state::facePickIndStart = mesh->nVertices();
            polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
            polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();
        }
        ImGui::SameLine(200);
        ImGui::DragFloat("Min Angle", &angleThresh, 0.05f, 0.05f, 0.5f, "%.2f");

        if (subSection("RemeshingOptions"))
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

            // Read from Options each time, so a switch of run shows the options of the new one
            float sizeMin = Options.min_absolute_length;
            float sizeMax = Options.max_absolute_length;
            float aspectMin = Options.aspect_min;
            float refineAngle = Options.refine_angle;
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
        if (subSection("Inspection"))
        {

            if (ImGui::Button("Bending Energy"))
            {
                psMesh->addVertexScalarQuantity("Bending Energy", M3DG.Sim_handler->Ev_Bending(std::vector<double>{1.0, 0.0}));
            }
            if (ImGui::Button("Bending Energy tan"))
            {
                psMesh->addVertexScalarQuantity("Bending Energy Tan", M3DG.Sim_handler->Ev_Bending_tan(std::vector<double>{1.0, 0.0}));
            }

            ImGui::TreePop();
        }

        ImGui::TreePop(); // This is required at the end of the if block
    }

    if (menuSection("Geometry"))
    {
        if (ImGui::Button("Vertex valence"))
        {
            VertexData<int> Valence(*mesh, 0);
            for (Vertex v : mesh->vertices())
            {
                for (Vertex v1 : v.adjacentVertices())
                    Valence[v] += 1;
            }

            psMesh->addVertexScalarQuantity("Valence", Valence);
        }

        if (ImGui::Button("Mean Curvature"))
        {
            VertexData<double> H(*mesh, 0.0);
            for (Vertex v : mesh->vertices())
                H[v] = geometry->scalarMeanCurvature(v);
            psMesh->addVertexScalarQuantity("H", H);
        }

        // Principal curvatures k_max >= k_min of the current shape (library values,
        // 1/r on a sphere of radius r). The field goes on the mesh and the extreme
        // value next to the button.
        static double max_curvature = 0.0, min_curvature = 0.0;
        static bool have_max_curvature = false, have_min_curvature = false;
        if (ImGui::Button("Max Curvature"))
        {
            geometry->refreshQuantities(); // the positions moved since the last request
            geometry->requireVertexMaxPrincipalCurvatures();
            VertexData<double> k_max = geometry->vertexMaxPrincipalCurvatures;
            geometry->unrequireVertexMaxPrincipalCurvatures();
            max_curvature = -std::numeric_limits<double>::infinity();
            for (Vertex v : mesh->vertices())
                max_curvature = std::max(max_curvature, k_max[v]);
            have_max_curvature = true;
            psMesh->addVertexScalarQuantity("Max curvature", k_max);
        }
        if (have_max_curvature)
        {
            ImGui::SameLine();
            ImGui::Text("max = %.4f", max_curvature);
        }
        if (ImGui::Button("Min Curvature"))
        {
            geometry->refreshQuantities();
            geometry->requireVertexMinPrincipalCurvatures();
            VertexData<double> k_min = geometry->vertexMinPrincipalCurvatures;
            geometry->unrequireVertexMinPrincipalCurvatures();
            min_curvature = std::numeric_limits<double>::infinity();
            for (Vertex v : mesh->vertices())
                min_curvature = std::min(min_curvature, k_min[v]);
            have_min_curvature = true;
            psMesh->addVertexScalarQuantity("Min curvature", k_min);
        }
        if (have_min_curvature)
        {
            ImGui::SameLine();
            ImGui::Text("min = %.4f", min_curvature);
        }

        ImGui::TreePop();
    }

    ImGuiTreeNodeFlags flag3 = ImGuiTreeNodeFlags_None;
    // Ok i want to do the energies
    // so the thing is that every energy has a different function that gets called. And i want the correct one to be called.
    if (menuSection("Energy"))
    {
        if (ImGui::Button("Display energy"))
        {
            double E = 0;
            Sim_handler.Calculate_energies(&E);
            std::cout << "THe current energy is  " << E << " \n";
        }
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

        static std::vector<std::vector<float>> Menu_E_constants;
        Menu_E_constants.resize(Energy_constants.size()); // another run can have other terms
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
            if (subSection(E_name.c_str()))
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
            if (subSection(("Bead " + std::to_string(bi)).c_str()))
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
    if (menuSection("Transformation"))
    {
        if (subSection("Rescale"))
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

    if (menuSection("Selection"))
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
            psMesh = registerMembrane();

            polyscope::state::facePickIndStart = mesh->nVertices();
            polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
            polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();
        }
        if (ImGui::Button("Flip Edge selected"))
        {
            for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
                 it != polyscope::state::subset.edges.end(); ++it)
            {
                Edge e = mesh->edge(*it);
                // Ok here i want to flip the edge
                MutationManager mm(*mesh, *geometry);
                mm.flipEdge(e);
            }
            psMesh = registerMembrane();
        }
        if (ImGui::Button("Collapse Edge selected"))
        {
            for (std::set<size_t>::iterator it = polyscope::state::subset.edges.begin();
                 it != polyscope::state::subset.edges.end(); ++it)
            {
                std::cout << "Collapsing one edge\n";
                Edge e = mesh->edge(*it);
                // Ok here i want to flip the edge
                MutationManager mm(*mesh, *geometry);
                mm.collapseEdge(e, 0.5);

                geometry->refreshQuantities();
                mesh->compress();
                geometry->inputVertexPositions = geometry->vertexPositions;
            }
            psMesh = registerMembrane();
        }
        if (ImGui::Button("SaveMesh"))
        {
            // I want to save the mesh.
            Save_mesh(basic_name, current_t);
        }

        if (ImGui::Button("FlipDelaunay"))
        {

            MutationManager mm(*mesh, *geometry);
            fixDelaunay(*mesh, *geometry, mm, Options);
        }
        ImGui::TreePop();
    }
    if (menuSection("Dump"))
    {

        if (ImGui::Button("Edge data"))
        {
            // OK so we dump the data
            std::ofstream f(basic_name + "Edge_data_step_" + std::to_string(current_t) + ".txt", std::ios_base::app);
            double H_e;
            for (Edge e : mesh->edges())
            {
                // We need to calculate some of this stuff
                H_e = geometry->scalarMeanCurvature(e.halfedge().vertex()) + geometry->scalarMeanCurvature(e.halfedge().twin().vertex());
                H_e /= 2;
                f << e.getIndex() << " " << geometry->edgeLength(e) << " " << geometry->dihedralAngle(e.halfedge()) << " " << H_e << "\n";
            }
            f.close();
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

            psMesh = registerMembrane();
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
        psMesh = registerMembrane();
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
        psMesh = registerMembrane();

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

bool ask_yes_no(const std::string &question)
{
    std::cout << question << " [y/N] " << std::flush;
    std::string answer;
    if (!std::getline(std::cin, answer))
        return false;
    return !answer.empty() && (answer[0] == 'y' || answer[0] == 'Y');
}

// New run folder holding the current state, without stepping: Input_file.json,
// membrane_0.obj and one row (step 0) in Output_data.txt and each bead file,
// so the folder can be browsed and continued like one written by main_cluster.
void create_run_folder(const std::string &config_path, const std::string &dir)
{
    make_dirs(dir);
    {
        std::ifstream src(config_path, std::ios::binary);
        std::ofstream dst(dir + "Input_file.json", std::ios::binary);
        dst << src.rdbuf();
    }
    record_run(sim.cfg.first_dir, dir, build_descriptive_name(sim), config_path);

    std::vector<std::string> bead_files = open_output_files(sim, dir, false);
    Save_mesh(mesh, geometry, dir, 0);

    Sim_handler.Calculate_gradient();
    double tot_E;
    Sim_handler.Calculate_energies(&tot_E);
    geometry->requireFaceAreas();
    double A = geometry->totalArea();
    geometry->unrequireFaceAreas();
    double V = geometry->totalVolume();

    std::ofstream out(dir + "Output_data.txt", std::ios_base::app);
    M3DG.write_output_row(out, 0.0, V, A, tot_E, 0.0);
    M3DG.write_bead_rows(bead_files);
    std::cout << "Created " << dir << " with the initial state, E = " << tot_E << "\n";
}

int main(int argc, char **argv)
{
    // --load_chain: the folder is a chain, its link*/<n> run folders are the runs of the menu. The rest
    // is the command line of the other programs
    bool load_chain = false;
    std::vector<char *> kept_args;
    for (int i = 0; i < argc; i++)
        if (i > 0 && std::string(argv[i]) == "--load_chain")
            load_chain = true;
        else
            kept_args.push_back(argv[i]);
    argc = int(kept_args.size());
    kept_args.push_back(nullptr);
    argv = kept_args.data();

    std::string config_path;
    bool check_only;
    if (!parse_command_line(argc, argv, config_path, check_only))
    {
        std::cerr << "       " << argv[0] << " --load_chain <chain folder: link*/<n> run folders>\n";
        return EXIT_FAILURE;
    }

    try
    {
        // A directory instead of an input file: its sub-folders with an
        // Input_file.json are the runs of the menu, and the first one opens
        // (with --load_chain the ones inside its link folders)
        struct stat path_info;
        const bool is_dir = stat(config_path.c_str(), &path_info) == 0 && S_ISDIR(path_info.st_mode);
        if (load_chain && !is_dir)
            throw std::runtime_error("--load_chain needs the folder of the chain, " + config_path + " is not a folder");
        SimConfig cfg;
        if (is_dir)
        {
            std::string root = config_path;
            if (root.back() != '/')
                root += '/';
            if (load_chain)
                chain_root = root;
            scan_run_folders_in(root, load_chain);
            if (run_folders.empty())
            {
                if (load_chain)
                    throw std::runtime_error("no link folder of " + root + " has a run folder with an Input_file.json");
                // The runs are one level further down: a chain
                scan_run_folders_in(root, true);
                throw std::runtime_error("no sub-folder of " + root + " has an Input_file.json" +
                                         (run_folders.empty() ? "" : " (its link folders hold " + std::to_string(run_folders.size()) + " runs: --load_chain opens a chain)"));
            }
            std::cout << "Found " << run_folders.size() << " runs in " << root << ", opening " << run_folders[0].name << "\n";
            cfg = run_config(root + run_folders[0].name + "/");
        }
        else
            // With "Subfolder" the run in first_dir/Subfolder is reopened and continued
            cfg = load_config(config_path, true);
        if (cfg.subfolder_missing &&
            !ask_yes_no("Create " + cfg.subfolder_dir + " with the initial state of " + config_path + "?"))
            return EXIT_FAILURE;
        build_simulation(cfg, sim);
        if (cfg.subfolder_missing)
            create_run_folder(config_path, cfg.subfolder_dir);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error setting up the simulation: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    try
    {
        prepare_run();
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
    scan_run_folders();

    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

    polyscope::TransparencyMode::Pretty;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None; // no floor, for the figures
    polyscope::init();
    applyMenuStyle();
    applyBackground();

    std::cout << "Im here happy\n";

    // std::sort(mesh_names.begin(), mesh_names.end());

    // return 1;

    // polyscope::state::userCallback = Callback_curvatures;
    // polyscope::state::userCallback = Callback_anim;

    polyscope::state::userCallback = Callback_qts;

    geometry->requireVertexPositions();

    psMesh = registerMembrane();
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
        std::vector<Vector3> beadPositions;
        for (size_t i = 0; i < Beads.size(); i++)
            beadPositions.push_back(Beads[i].Pos);
        registerBeadCloud(beadPositions);
    }

    // std::cout<<"The radius would be " <<Beads[0].sigma*1.0 <<"\n";
    // psCloud->setPointRadius(Beads[0].sigma);
    // }

    std::cout << "Here7\n";


    // i JUST WANT SOMETHING ELSE HERE
    double lengthScale = geometry->meanEdgeLength();
    vertexRadius = 0.002;
    edgeRadius = 0.001;

    frameCameraOnMembrane();

    std::cout << "Calling the mesh\n";
    polyscope::show();

    return EXIT_SUCCESS;
}
