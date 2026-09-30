// Legacy Mem3DG code, not compiled.
//
// Moved out of src/Mem-3dg.cpp during the cleanup: the old per-energy force API,
// the integrators and line searches that took the energy constants as arguments,
// finite-difference checks nothing calls, and a copy of the remesher (the one in
// use lives in deps/geometry-central/src/surface/remeshing.cpp).
// Kept for reference; see git history for the original context.

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
VertexData<Vector3> Mem3DG::buildFlowOperator(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd)
{

  // Lets get our target area and curvature

  double V = geometry->totalVolume();
  double D_P = -P0 * (V - V_bar) / (V_bar * V_bar);

  double A_bar = 4 * PI * pow(3 * V_bar / (4 * PI * nu), 2.0 / 3.0);
  double H_bar = sqrt(4 * PI / A_bar) * c0 / 2.0; // Coment this with another comment
  double A = geometry->totalArea();
  double lambda = KA * (A - A_bar) / (A_bar * A_bar);

  // This lines are for the bunny i need to delete them later
  // lambda=KA;
  // KB=0;
  // return (SurfaceTension(lambda)+OsmoticPressure(D_P));

  return (KB * Bending(H_bar) + D_P * OsmoticPressure() + lambda * SurfaceTension());
  // return (Bending(H_bar,KB)+SurfaceTension(lambda));

  //
  // +SurfaceTension(lambda)
}

VertexData<Vector3> Mem3DG::buildFlowOperator(double V_bar, double P0, double KA, double KB, double h)
{

  // Lets get our target area and curvature

  double V = geometry->totalVolume();
  double D_P = -P0 * (V - V_bar) / V_bar / V_bar;

  // double A_bar=4*PI*pow(3*V_bar/(4*PI*nu),2.0/3.0);
  double H_bar = 0.0; // Coment this with another comment
  double A = geometry->totalArea();
  double lambda = KA;

  // This lines are for the bunny i need to delete them later
  // lambda=KA;
  // KB=0;
  // return (SurfaceTension(lambda)+OsmoticPressure(D_P));
  // return (Bending(H_bar,KB)+lambda*SurfaceTension());

  return (KB * Bending(H_bar) + D_P * OsmoticPressure() + lambda * SurfaceTension());
  // return (Bending(H_bar,KB)+SurfaceTension(lambda));

  //
  // +SurfaceTension(lambda)
}

VertexData<Vector3> Mem3DG::buildFlowOperator(double h, double V_bar, double P0, double KA)
{

  // Lets get our target area and curvature

  double V = geometry->totalVolume();
  double D_P = -P0 * (V - V_bar) / V_bar / V_bar;

  return (D_P * OsmoticPressure() + KA * SurfaceTension());
  // return (Bending(H_bar,KB)+SurfaceTension(lambda));

  //
  // +SurfaceTension(lambda)
}

VertexData<Vector3> Mem3DG::OsmoticPressure() const
{

  // You have the face normals

  size_t index;
  Vector3 Normal;
  size_t N_vert = mesh->nVertices();

  VertexData<Vector3> Force(*mesh);
  for (Vertex v : mesh->vertices())
  {
    // do science here
    index = v.getIndex();
    Normal = {0, 0, 0};
    for (Face f : v.adjacentFaces())
    {
      Normal += geometry->faceArea(f) * geometry->faceNormal(f);
    }
    // Force[v.getIndex()]=D_P*Normal/3.0;
    Force[v.getIndex()] = Normal / 3.0;
  }

  // std::cout<< "THe osmotic pressure force in magnitude is: "<< D_P*sqrt(Force.transpose()*Force) <<"\n";
  return Force;
}

VertexData<Vector3> Mem3DG::SurfaceTension() const
{

  size_t index;
  Vector3 Normal;
  size_t N_vert = mesh->nVertices();
  VertexData<Vector3> Force(*mesh);

  for (Vertex v : mesh->vertices())
  {
    Normal = {0, 0, 0};
    // for(Halfedge he: v.outgoingHalfedges()){
    //   Normal+=2*computeHalfedgeMeanCurvatureVector(he);
    // }
    // index=v.getIndex();

    Normal = 2 * geometry->vertexNormalMeanCurvature(v);
    Force[v] = Normal;
  }

  // std::cout<< "THe surface tension force in magnitude is: "<< -1*lambda*sqrt(Force.transpose()*Force) <<"\n";
  return -1 * Force;
}

Vector3 Mem3DG::computeHalfedgeMeanCurvatureVector(Halfedge he) const
{
  size_t fID = he.face().getIndex();
  size_t fID_he_twin = he.twin().face().getIndex();
  Vector3 areaGrad{0, 0, 0};

  Vector3 EdgeVector = geometry->inputVertexPositions[he.next().next().vertex()] - geometry->inputVertexPositions[he.next().vertex()];
  Vector3 EdgeVector2 = geometry->inputVertexPositions[he.twin().vertex()] - geometry->inputVertexPositions[he.twin().next().next().vertex()];

  areaGrad +=
      0.25 * cross(geometry->faceNormal(he.face()), EdgeVector);

  areaGrad += 0.25 * cross(geometry->faceNormal(he.twin().face()),
                           EdgeVector2);

  return areaGrad / 2;
}

Vector3 Mem3DG::computeHalfedgeGaussianCurvatureVector(Halfedge he) const
{
  Vector3 gaussVec{0, 0, 0};
  if (!he.edge().isBoundary())
  {
    // gc::Vector3 eji{} = -vecFromHalfedge(he, *vpg);
    gaussVec = 0.5 * geometry->dihedralAngle(he) * (-1 * geometry->inputVertexPositions[he.next().vertex()] + geometry->inputVertexPositions[he.vertex()]).unit();
  }
  else
  {
    gaussVec = 0.5 * geometry->dihedralAngle(he) * (-1 * geometry->inputVertexPositions[he.next().vertex()] + geometry->inputVertexPositions[he.vertex()]).unit();
    // std::cout<< "Dihedral angle "<<0.5 * geometry->dihedralAngle(he)<<"\n";
    // std::cout<<" Unit vector of an edge"<< ( -1* geometry->inputVertexPositions[he.next().vertex()]+geometry->inputVertexPositions[he.vertex()] ).unit()<<"\n";
    // std::cout<<"This mean gaussian curvature shouldnt work";
  }
  return gaussVec;
}

Vector3 Mem3DG::dihedralAngleGradient(Halfedge he, Vertex v) const
{
  // std::cout<< he.edge().isBoundary();

  double l = geometry->edgeLength(he.edge());

  if (he.edge().isBoundary())
  {
    return Vector3{0, 0, 0};
  }
  else if (he.vertex() == v)
  { // This is only used for the SIJ_1
    return (geometry->cotan(he.next().next()) *
                geometry->faceNormal(he.face()) +
            geometry->cotan(he.twin().next()) *
                geometry->faceNormal(he.twin().face())) /
           l;
  }
  else if (he.next().vertex() == v)
  { // This is for the firt s term
    return (geometry->cotan(he.twin().next().next()) *
                geometry->faceNormal(he.twin().face()) +
            geometry->cotan(he.next()) *
                geometry->faceNormal(he.face())) /
           l;
  }
  else if (he.next().next().vertex() == v)
  { // Este ocurre para el segundo termino
    return (-(geometry->cotan(he.next().next()) +
              geometry->cotan(he.next())) *
            geometry->faceNormal(he.face())) /
           l;
  }
  else
  {
    // mem3dg_runtime_error("Unexpected combination of halfedge and vertex!");
    std::cout << "THe dihedral angle gradient is not working\n";
    return Vector3{0, 0, 0};
  }
  std::cout << "THis is impossible to print\n";
  return Vector3{0, 0, 0};
}

VertexData<Vector3> Mem3DG::Bending(double H0) const
{

  size_t neigh_index;
  size_t N_vert = mesh->nVertices();

  Vector3 Hij;
  Vector3 Kij;
  Vector3 Sij_1;
  Vector3 Sij_2;

  Vector3 F1 = {0, 0, 0};
  Vector3 F2 = {0, 0, 0};
  Vector3 F3 = {0, 0, 0};
  Vector3 F4 = {0, 0, 0};

  Vector3 Position_1;
  Vector3 Position_2;

  VertexData<Vector3> Force(*mesh);

  VertexData<double> Scalar_MC(*mesh, 0.0);
  double factor;
  size_t index1;
  for (Vertex v1 : mesh->vertices())
  {
    index1 = v1.getIndex();
    Scalar_MC[index1] = geometry->scalarMeanCurvature(v1) / geometry->barycentricDualArea(v1);
  }

  auto start = chrono::steady_clock::now();
  double H0i;
  double H0j;
  size_t index;
  for (Vertex v : mesh->vertices())
  {
    F1 = {0, 0, 0};
    F2 = {0, 0, 0};
    F3 = {0, 0, 0};
    F4 = {0, 0, 0};
    index = v.getIndex();
    // H0i= (system_time<50? (H_Vector_0[index]+H0)/2.0: H0);
    H0i = H0;
    Position_1 = geometry->inputVertexPositions[v];
    for (Halfedge he : v.outgoingHalfedges())
    {

      neigh_index = he.tipVertex().getIndex();
      // H0j= (system_time<50? (H_Vector_0[neigh_index]+H0)/2: H0);
      H0j = H0;

      Position_2 = geometry->inputVertexPositions[neigh_index];

      Kij = computeHalfedgeGaussianCurvatureVector(he);
      // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))-1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
      factor = -1 * (Scalar_MC[index] - H0i) - 1 * (Scalar_MC[neigh_index] - H0j);

      F1 = F1 + factor * Kij;

      Hij = 2 * computeHalfedgeMeanCurvatureVector(he);

      // factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0)*(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0)*(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
      factor = (1 / 3.0) * (Scalar_MC[index] * Scalar_MC[index] - H0i * H0i) + (2.0 / 3.0) * (Scalar_MC[neigh_index] * Scalar_MC[neigh_index] - H0j * H0j);

      F2 = F2 + factor * Hij;

      Sij_1 = geometry->edgeLength(he.edge()) * dihedralAngleGradient(he, he.vertex());

      // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0));
      factor = -1 * (Scalar_MC[index] - H0i);

      F3 = F3 + factor * Sij_1;

      Sij_2 = (geometry->edgeLength(he.twin().edge()) * dihedralAngleGradient(he.twin(), he.vertex()) + geometry->edgeLength(he.next().edge()) * dihedralAngleGradient(he.next(), he.vertex()) + geometry->edgeLength(he.twin().next().next().edge()) * dihedralAngleGradient(he.twin().next().next(), he.vertex()));

      // Sij_2=-1*( geometry->cotan(he.next().next())*geometry->faceNormal(he.face()) + geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()));

      // factor= -1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
      factor = -1 * (Scalar_MC[neigh_index] - H0j);

      F4 = F4 + factor * Sij_2;
    }

    Force[index] = F1 + F2 + F3 + F4;
  }

  return Force;
}

//
double Mem3DG::E_Volume_constraint(double KV, double V, double V_bar) const
{
  return 0.5 * KV * (V - V_bar) * (V - V_bar) / (V_bar * V_bar);
}

double Mem3DG::E_Pressure(double P0, double V, double V_bar) const
{

  // double V = geometry->totalVolume();

  return -1 * 0.5 * P0 * (V - V_bar);
}

double Mem3DG::E_Area_constraint(double KA, double A, double A_bar) const
{

  return 0.5 * KA * (A - A_bar) * (A - A_bar) / (A_bar * A_bar);
}

double Mem3DG::E_Surface(double KA, double A, double A_bar) const
{

  // return 0.5*KA*A*A;
  return 0.5 * KA * (A - A_bar) * (A - A_bar) / (A_bar * A_bar);
}

double Mem3DG::E_Bending(double H0, double KB) const
{
  size_t index;
  double Eb = 0;
  double H;
  double r_eff2;
  Vector3 Pos;
  for (Vertex v : mesh->vertices())
  {
    // boundary_fix
    if (v.isBoundary())
      continue;
    index = v.getIndex();
    // Scalar_MC.coeffRef(index)
    Pos = geometry->inputVertexPositions[v];
    r_eff2 = Pos.z * Pos.z + Pos.y * Pos.y;
    if (r_eff2 > 1.6 && boundary)
      continue;

    H = abs(geometry->scalarMeanCurvature(v) / geometry->barycentricDualArea(v));

    if (std::isnan(H))
    {
      continue;
      std::cout << "Dual area: " << geometry->barycentricDualArea(v);
      std::cout << "Scalar mean Curv" << geometry->scalarMeanCurvature(v);
      std::cout << "One of the H is not a number\n";
    }

    Eb += KB * H * H * geometry->barycentricDualArea(v);
  }

  return Eb;
}

SparseMatrix<double> Mem3DG::H2_operator(bool CM, bool Vol_const, bool Area_const)
{

  // Ok so i have the gradient i want to calculate

  // std::cout<<"We are doing the sobolev operator ! \n";
  // std::cout<<"THe constraints are CM: "<< CM << " Volume " <<  Vol_const <<" and area " << Area_const <<" \n";
  SparseMatrix<double> L = geometry->laplaceMatrix();
  SparseMatrix<double> M = geometry->massMatrix();
  SparseMatrix<double> Inv_M(mesh->nVertices(), mesh->nVertices());
  for (size_t index = 0; index < mesh->nVertices(); index++)
  {
    Inv_M.coeffRef(index, index) = 1.0 / (M.coeff(index, index));
  }

  SparseMatrix<double> J = L.transpose() * Inv_M * L;
  // I need the constraint gradientsx
  VertexData<Vector3> grad_sur = SurfaceGrad();
  VertexData<Vector3> grad_vol = OsmoticPressure();
  // VertexData<Vector3> grad_ben = Bending(0.0);

  // i CAN MAYBE MULTIPLY BY THE SURFACE AREAS

  // std::cout<<"Grad ben \n";

  // std::cout<<"We are debugginggg \n";
  size_t N_vert = mesh->nVertices();

  int Num_constraints = 0;
  if (CM)
    Num_constraints += 3;
  if (Vol_const)
    Num_constraints += 1;
  if (Area_const)
    Num_constraints += 1;

  // std::cout<<"THe number of constraints is "<< Num_constraints<<" \n";

  SparseMatrix<double> S(N_vert * 3 + Num_constraints, N_vert * 3 + Num_constraints);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  int highest_row = 0;
  int highest_col = 0;
  // We add the constraints first

  int Constraint_number = 0;

  for (size_t index = 0; index < mesh->nVertices(); index++)
  {

    Constraint_number = 0;
    //  Area constraint
    if (Area_const)
    {
      // std::cout<<"Area constraint on \n";
      tripletList.push_back(T(3 * index, 3 * N_vert + Constraint_number, grad_sur[index].x));
      tripletList.push_back(T(3 * index + 1, 3 * N_vert + Constraint_number, grad_sur[index].y));
      tripletList.push_back(T(3 * index + 2, 3 * N_vert + Constraint_number, grad_sur[index].z));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index, grad_sur[index].x));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 1, grad_sur[index].y));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 2, grad_sur[index].z));

      Constraint_number++;
    }

    // Volume Constraint
    if (Vol_const)
    {
      // std::cout<<"Vol const on\n";
      tripletList.push_back(T(3 * index, 3 * N_vert + Constraint_number, grad_vol[index].x));
      tripletList.push_back(T(3 * index + 1, 3 * N_vert + Constraint_number, grad_vol[index].y));
      tripletList.push_back(T(3 * index + 2, 3 * N_vert + Constraint_number, grad_vol[index].z));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index, grad_vol[index].x));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 1, grad_vol[index].y));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 2, grad_vol[index].z));

      Constraint_number++;
    }
    // Position Constraint I

    if (CM)
    {
      // std::cout<<"POs constraint on \n";
      tripletList.push_back(T(3 * index, 3 * N_vert + Constraint_number, 1));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index, 1));
      Constraint_number++;

      tripletList.push_back(T(3 * index + 1, 3 * N_vert + Constraint_number, 1));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 1, 1));
      Constraint_number++;

      tripletList.push_back(T(3 * index + 2, 3 * N_vert + Constraint_number, 1));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 2, 1));
      Constraint_number++;
    }
  }

  //  std::cout<<"Before setting from tripplets\n";
  // std::cout<<"THe number of vertices is " << N_vert <<"\n";
  // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
  // std::cout<<"THe calculated col is " << highest_col <<" \n";
  // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
  // std::cout<<"THe calculated row is " << highest_row <<" \n";
  // Now we iterate over the Laplacian
  int row;
  int col;
  double value;

  for (long int k = 0; k < J.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(J, k); it; ++it)
    {
      value = it.value();
      row = it.row();
      col = it.col();
      tripletList.push_back(T(3 * row, 3 * col, value));
      tripletList.push_back(T(3 * row + 1, 3 * col + 1, value));
      tripletList.push_back(T(3 * row + 2, 3 * col + 2, value));
      // if( 3*row > highest_row) highest_row = 3*row;
      // if( 3*col > highest_col) highest_col = 3*col;
    }
  }

  // std::cout<<"Before setting from tripplets\n";
  // std::cout<<"THe number of vertices is " << N_vert <<"\n";
  // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
  // std::cout<<"THe calculated col is " << highest_col <<" \n";
  // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
  // std::cout<<"THe calculated row is " << highest_row <<" \n";
  S.setFromTriplets(tripletList.begin(), tripletList.end());
  // std::cout<<"THe matrix has been formed\n";

  return S;

  // std::cout<<"AFTER setting from tripplets\n";
  // // We need to create the vector of the RHS

  // Vector<double> RHS(N_vert*3+Num_constraints);
  // int highest_idx = 0;
  // for(size_t index; index<mesh->nVertices();index++)
  // {
  //     RHS.coeffRef(3*index)=grad_ben[index].x;
  //     RHS.coeffRef(3*index+1)=grad_ben[index].y;
  //     RHS.coeffRef(3*index+2)=grad_ben[index].z;
  //     if( 3*index +2 > highest_idx) highest_idx = 3*index+2;
  // }

  // // std::cout<<"The highest index is " << highest_idx <<"\n";
  // // std::cout<<"CREATING RHS\n";
  // for(size_t index = 0; index<Num_constraints;index++)
  // {
  //     RHS.coeffRef(3*N_vert+index)=0;
  // }
  // // RHS.coeffRef(3*N_vert)=0;
  // // RHS.coeffRef(3*N_vert+1)=0;
  // // RHS.coeffRef(3*N_vert+2)=0;
  // // RHS.coeffRef(3*N_vert+3)=0;
  // // RHS.coeffRef(3*N_vert+4)=0;

  // // std::cout<<"RHS IMPLEMENTED\n";

  // // std::cout<<"The RHS is " << RHS << "\n";
  // // We have The RHS and the matrix, its tiiiime
  // // std::cout<<"Lets solve \n";
  // Eigen::SparseLU<SparseMatrix<double>> solver;
  // // std::cout<<"Solver defined \n";

  // solver.compute(S);

  // Vector<double> result = solver.solve(RHS);
  // // std::cout<<"Result retrieved\n";
  // VertexData<Vector3> Final_Force(*mesh);
  // for(size_t index; index<mesh->nVertices();index++)
  // {
  //     Final_Force[index]=Vector3{result.coeff(3*index),result.coeff(3*index+1),result.coeff(3*index+2)};
  // }

  // // std::cout<<"The two lambda for the constraints are "<< result[N_vert] << " and "<< result[N_vert+1] <<"\n";
  // // return Final_Force;
}

SparseMatrix<double> Mem3DG::H1_operator(bool CM, bool Vol_const, bool Area_const)
{

  // Ok so i have the gradient i want to calculate

  // std::cout<<"We are doing the sobolev operator ! \n";
  // std::cout<<"THe constraints are CM: "<< CM << " Volume " <<  Vol_const <<" and area " << Area_const <<" \n";
  SparseMatrix<double> L = geometry->laplaceMatrix();
  // SparseMatrix<double> M = geometry->massMatrix();
  // SparseMatrix<double> Inv_M(mesh->nVertices(),mesh->nVertices());
  // for(size_t index= 0; index<mesh->nVertices();index++)
  // {
  //     Inv_M.coeffRef(index,index)= 1.0/(M.coeff(index,index));

  // }

  SparseMatrix<double> J = L;
  // I need the constraint gradientsx
  VertexData<Vector3> grad_sur = SurfaceGrad();
  VertexData<Vector3> grad_vol = OsmoticPressure();
  // VertexData<Vector3> grad_ben = Bending(0.0);

  // i CAN MAYBE MULTIPLY BY THE SURFACE AREAS

  // std::cout<<"Grad ben \n";

  // std::cout<<"We are debugginggg \n";
  size_t N_vert = mesh->nVertices();

  int Num_constraints = 0;
  if (CM)
    Num_constraints += 3;
  if (Vol_const)
    Num_constraints += 1;
  if (Area_const)
    Num_constraints += 1;

  // std::cout<<"THe number of constraints is "<< Num_constraints<<" \n";

  SparseMatrix<double> S(N_vert * 3 + Num_constraints, N_vert * 3 + Num_constraints);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  int highest_row = 0;
  int highest_col = 0;
  // We add the constraints first

  int Constraint_number = 0;

  for (size_t index = 0; index < mesh->nVertices(); index++)
  {

    Constraint_number = 0;
    //  Area constraint
    if (Area_const)
    {
      // std::cout<<"Area constraint on \n";
      tripletList.push_back(T(3 * index, 3 * N_vert + Constraint_number, grad_sur[index].x));
      tripletList.push_back(T(3 * index + 1, 3 * N_vert + Constraint_number, grad_sur[index].y));
      tripletList.push_back(T(3 * index + 2, 3 * N_vert + Constraint_number, grad_sur[index].z));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index, grad_sur[index].x));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 1, grad_sur[index].y));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 2, grad_sur[index].z));

      Constraint_number++;
    }

    // Volume Constraint
    if (Vol_const)
    {
      // std::cout<<"Vol const on\n";
      tripletList.push_back(T(3 * index, 3 * N_vert + Constraint_number, grad_vol[index].x));
      tripletList.push_back(T(3 * index + 1, 3 * N_vert + Constraint_number, grad_vol[index].y));
      tripletList.push_back(T(3 * index + 2, 3 * N_vert + Constraint_number, grad_vol[index].z));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index, grad_vol[index].x));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 1, grad_vol[index].y));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 2, grad_vol[index].z));

      Constraint_number++;
    }
    // Position Constraint I

    if (CM)
    {
      // std::cout<<"POs constraint on \n";
      tripletList.push_back(T(3 * index, 3 * N_vert + Constraint_number, 1));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index, 1));
      Constraint_number++;

      tripletList.push_back(T(3 * index + 1, 3 * N_vert + Constraint_number, 1));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 1, 1));
      Constraint_number++;

      tripletList.push_back(T(3 * index + 2, 3 * N_vert + Constraint_number, 1));
      tripletList.push_back(T(3 * N_vert + Constraint_number, 3 * index + 2, 1));
      Constraint_number++;
    }
  }

  //  std::cout<<"Before setting from tripplets\n";
  // std::cout<<"THe number of vertices is " << N_vert <<"\n";
  // std::cout<<"THe highest expected col is " << 3*N_vert+2 <<" \n";
  // std::cout<<"THe calculated col is " << highest_col <<" \n";
  // std::cout<<"THe highest expected row is " << 3*N_vert+2 <<" \n";
  // std::cout<<"THe calculated row is " << highest_row <<" \n";
  // Now we iterate over the Laplacian
  int row;
  int col;
  double value;

  for (long int k = 0; k < J.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(J, k); it; ++it)
    {
      value = it.value();
      row = it.row();
      col = it.col();
      tripletList.push_back(T(3 * row, 3 * col, value));
      tripletList.push_back(T(3 * row + 1, 3 * col + 1, value));
      tripletList.push_back(T(3 * row + 2, 3 * col + 2, value));
      // if( 3*row > highest_row) highest_row = 3*row;
      // if( 3*col > highest_col) highest_col = 3*col;
    }
  }

  S.setFromTriplets(tripletList.begin(), tripletList.end());

  return S;
}

VertexData<Vector3> Mem3DG::Linear_force_field(double x0, double slope) const
{
  // Can i make it volume preserving?

  VertexData<Vector3> FField(*mesh, Vector3({0, 0, 0}));
  VertexData<Vector3> FField_correction(*mesh, Vector3({0, 0, 0}));
  Vector3 Normal;
  double mag = 0;
  double tot_mag = 0;
  for (Vertex v : mesh->vertices())
  {
    // I want to iterate over all vertices
    Normal = geometry->vertexNormalAreaWeighted(v);
    mag = abs(geometry->inputVertexPositions[v].x - x0) * slope;
    FField[v] = mag * Normal;
    tot_mag += mag;
    FField_correction[v] = Normal;
  }

  // int index;
  // for(Vertex v : mesh->vertices()) {
  //    // do science here
  //       index=v.getIndex();
  //       Normal={0,0,0};
  //       for(Face f : v.adjacentFaces()) {
  //           Normal+=geometry->faceArea(f)*geometry->faceNormal(f);

  //       }
  //       // Force[v.getIndex()]=D_P*Normal/3.0;
  //       mag = abs(geometry->inputVertexPositions[v].x-x0)*slope;
  //       tot_mag += mag;
  //       FField[v.getIndex()]=mag*Normal/3.0;
  //   }

  tot_mag = -1 * tot_mag / mesh->nVertices();

  // Now that i have the area average of this magnitude we can add it
  for (Vertex v : mesh->vertices())
  {
    FField_correction[v] *= tot_mag;
  }

  // FField_correction = OsmoticPressure(tot_mag);
  // return FField;
  return FField + FField_correction;
}

void Mem3DG::Smooth_vertices()
{

  size_t N_vert = mesh->nVertices();
  SparseMatrix<double> L(N_vert, N_vert);

  L = geometry->uniformlaplacianMatrix();

  // Ok i have the matrix now what
  Vector<double> xpos = Vector<double>(N_vert);
  Vector<double> ypos = Vector<double>(N_vert);
  Vector<double> zpos = Vector<double>(N_vert);
  Vector3 Pos;
  for (size_t i = 0; i < N_vert; i++)
  {
    Pos = geometry->inputVertexPositions[i];
    xpos[i] = Pos.x;
    ypos[i] = Pos.y;
    zpos[i] = Pos.z;
  }
  xpos = xpos + L * xpos;
  ypos = ypos + L * ypos;
  zpos = zpos + L * zpos;

  for (size_t i = 0; i < N_vert; i++)
  {
    geometry->inputVertexPositions[i] = Vector3({xpos[i], ypos[i], zpos[i]});
  }

  return;
}

double Mem3DG::Backtracking_grad_Normal_2(Eigen::VectorXd p_lambda, double Projection, double Current_grad_norm)
{

  double c1 = 0.0001;
  double rho = 0.5;
  double alpha = 1;
  // alpha = 5e-4;
  double position_Projeection = 0;
  double X_pos;

  int N_vert = mesh->nVertices();
  int N_beads = Beads.size();

  VertexData<Vector3> Step_newton = Sim_handler->Current_grad;
  std::vector<Vector3> Step_beads(0);

  for (size_t bi = 0; bi < N_beads; bi++)
    Step_beads.push_back(Beads[bi]->Total_force);

  double PrevNorm = Current_grad_norm;

  double NewNorm;
  double NewMerit;
  double V;
  double A;
  double PrevMerit = 0;
  Sim_handler->Calculate_Merit(&PrevMerit);

  // Now we have the whole merit function this can be turned into a proper function later

  VertexData<Vector3> initial_pos(*mesh);
  Eigen::VectorXd initial_lag;
  std::vector<Vector3> initial_bead_pos(0);

  // We will open a file to log the backtracking process
  std::ofstream backtrack_log;

  backtrack_log.open(basic_name + "backtrack_log.txt", std::ios::app);
  // So the order will be : Prevnorm Projection NEWNORM NEWNORM NEWNORM ...

  backtrack_log << PrevMerit << " ";
  if (recentering)
  {
    if (Field != "None")
    {
      // std::cout<<"THere is a field\n";
      Vector3 leftmost = Vector3({1e10, 1e10, 1e10});
      for (Vertex v : mesh->vertices())
      {
        if (geometry->inputVertexPositions[v].x < leftmost.x)
          leftmost = geometry->inputVertexPositions[v];
      }
      // I have the leftmost vertex, the position of this vertex should be P
      leftmost = Vector3({-1.0, 0.0, 0.0}) - leftmost;
      VertexData<Vector3> Displacement(*mesh, leftmost);
      geometry->inputVertexPositions += Displacement;
    }
    else
    {
      // Ok so if there is a field the way we normalize is different.
      //
      Vector3 CoM = geometry->centerOfMass();

      geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
    }
  }
  Vector3 CoM = geometry->centerOfMass();

  initial_pos = geometry->inputVertexPositions;
  initial_lag = Sim_handler->Lagrange_mult;

  for (int bi = 0; bi < N_beads; bi++)
  {
    initial_bead_pos.push_back(Beads[bi]->Pos);
  }

  std::vector<Vector3> Bead_init;

  for (int i = 0; i < N_beads; i++)
    Bead_init.push_back(Beads[i]->Pos);

  Vector3 center;
  geometry->inputVertexPositions += alpha * Step_newton;
  // We move the beads;
  for (size_t i = 0; i < Beads.size(); i++)
    Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));

  for (int i = 0; i < Sim_handler->N_constraints; i++)
  {
    if (Sim_handler->Constraints[i] == "Volume" || Sim_handler->Constraints[i] == "Area")
    {
      Sim_handler->Lagrange_mult[i] -= alpha * p_lambda[i];
    }
  }

  bool nanflag = false;
  for (Vertex v : mesh->vertices())
  {
    if (isnan(geometry->inputVertexPositions[v].norm2()))
      nanflag = true;
  }

  if (nanflag)
    std::cout << "At least one vertex has nan position\n";

  center = geometry->centerOfMass();
  Vector3 Vertex_pos;

  size_t bead_count = 0;
  NewNorm = 0.0;
  Sim_handler->Calculate_Lag_norm_Normal(&NewNorm);
  NewMerit = 0.0;
  Sim_handler->Calculate_Merit(&NewMerit);

  size_t counter = 0;
  bool displacement_cond = true;

  if (Projection < 1e-7)
  {
    small_TS = true;
    std::cout << "The norm  diff is quite small and so is the gradient\n";
    std::cout << "The norm diff is" << abs(NewNorm - PrevNorm) / PrevNorm << "\n";
    std::cout << "The projection is" << Projection << "\n";
    return -1.0;
  }

  // Ok so here we have the new merit and the

  while (true)
  {
    displacement_cond = true;
    backtrack_log << NewMerit << " ";

    for (size_t i = 0; i < Beads.size(); i++)
      displacement_cond = displacement_cond && Beads[i]->Total_force.norm() * alpha < 0.1 * Beads[i]->sigma;

    // if( fabs(PrevNorm-NewNorm) <= alpha * Projection && NewNorm < PrevNorm  ) {
    // if (NewNorm < PrevNorm)
    if (NewMerit < PrevMerit)
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
      // std::cout << "The merit diff is" << abs(NewMerit - PrevMerit) << "\n";

      // std::cout << "THe relative energy diff  is" << abs((NewNorm - PrevNorm) / PrevNorm) << "\n";
      // std::cout << "THe relative energy diff  is" << abs((NewMerit - PrevMerit) / PrevMerit) << "\n";
      // std::cout << " THe merit function is " << NewMerit << " \n";
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

    // geometry->refreshQuantities();

    bead_count = 0;

    // NewE = 0.0;
    // Sim_handler->Calculate_energies(&NewE);
    NewNorm = 0.0;
    Sim_handler->Calculate_Lag_norm_Normal(&NewNorm);
    NewMerit = 0.0;
    Sim_handler->Calculate_Merit(&NewMerit);
    // Sim_handler->Calculate_Lag_norm(&NewNorm);
    // NewNorm = NewNorm;
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
    {
      std::cout << "NotRecentering after crisis\n";
    }
    else
    {
      // std::cout<<"rENORMALIZING\n";
      CoM = geometry->centerOfMass();

      if (Field != "None")
      {
        // std::cout<<"THere is a field\n";
        Vector3 leftmost = Vector3({1e10, 1e10, 1e10});
        for (Vertex v : mesh->vertices())
        {
          if (geometry->inputVertexPositions[v].x < leftmost.x)
            leftmost = geometry->inputVertexPositions[v];
        }
        // I have the leftmost vertex, the position of this vertex should be P
        leftmost = Vector3({-1.0, 0.0, 0.0}) - leftmost;
        CoM = leftmost;
        VertexData<Vector3> Displacement(*mesh, leftmost);
        geometry->inputVertexPositions += Displacement;
      }
      else
      {
        geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
      }
      // CoM = geometry->centerOfMass();

      for (size_t i = 0; i < Beads.size(); i++)
      {

        // Here i need to move the bead
        // if(Beads[i]->state!= "froze"){
        Beads[i]->Pos -= CoM;
        // }
      }
      // }
      //
    }
    // geometry->refreshQuantities();
  }

  // std::cout<<"The difference in energy is " << fabs(NewE-previousE) <<"(: \n";
  // std::cout<<"The new norm is " << NewNorm << "\n";
  return alpha;
}

double Mem3DG::Backtracking(VertexData<Vector3> Force, double P0, double V_bar, double A_bar, double KA, double KB, double H_bar, bool bead, bool pulling)
{

  // std::cout<<"Backtracking\n";
  bool other_pulling = false;
  bool two_bead_pulling = true;
  double Bead_force_prev = -12;
  double Bead_force_new = -15;
  double c1 = 1e-4;
  double rho = 0.5;
  double alpha = 1e-3;
  double positionProjection = 0;
  double X_pos;

  // A=geometry->totalArea();
  // V=geometry->totalVolume();
  // E_Vol = E_Pressure(D_P,V,V_bar);
  // E_Sur = E_Surface(KA,A,A_bar);
  // E_Ben = E_Bending(H_bar,KB);
  // E_Bead = Bead_1.Energy();
  double previousE = E_Vol + E_Sur + E_Ben + E_Bead;
  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  if (recentering)
  {
    geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
  }
  initial_pos = geometry->inputVertexPositions;

  std::vector<Vector3> Bead_init;

  for (size_t i = 0; i < Beads.size(); i++)
    Bead_init.push_back(Beads[i]->Pos);

  // std::cout<<"THe current energy is "<<previousE <<"Is this awful?\n";
  // Zeroth iteration
  double Projection = 0;
  Vector3 center;

  geometry->inputVertexPositions += alpha * Force;
  // geometry->refreshQuantities();

  center = geometry->centerOfMass();
  Vector3 Vertex_pos;

  if (!pulling)
  {
    for (size_t i = 0; i < Beads.size(); i++)
      Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));
    // this->Bead_1.Move_bead(alpha,Vector3({0,0,0}));
  }

  Total_force = Vector3({0, 0, 0});
  for (Vertex v : mesh->vertices())
  {
    Projection += Force[v.getIndex()].norm2();
    Total_force += Force[v.getIndex()];
  }

  grad_norm = Projection;

  A = geometry->totalArea();
  V = geometry->totalVolume();
  double D_P = -1 * P0 * (V - V_bar) / (V_bar * V_bar);
  E_Vol = E_Pressure(D_P, V, V_bar);
  E_Sur = E_Surface(KA, A, A_bar);
  E_Ben = E_Bending(H_bar, KB);
  E_Bead = 0;
  for (size_t i = 0; i < Beads.size(); i++)
    E_Bead += Beads[i]->Energy();
  // E_Bead=Bead_1.Energy();
  NewE = E_Vol + E_Sur + E_Ben + E_Bead;

  if (std::isnan(E_Vol))
  {
    std::cout << "E vol is nan\n";
  }
  if (std::isnan(E_Sur))
  {
    std::cout << "E sur is nan\n";
  }
  if (std::isnan(E_Ben))
  {
    std::cout << "E ben is nan\n";
  }

  size_t counter = 0;

  bool displacement_cond = true;

  // std::cout<<"starting while\n";
  while (true)
  {
    // if(true){
    displacement_cond = true;
    for (size_t i = 0; i < Beads.size(); i++)
      displacement_cond = displacement_cond && Beads[i]->Total_force.norm() * alpha < 0.1;
    // std::cout<<"Displacement cond is"<< displacement_cond <<"\n";
    if (NewE <= previousE - c1 * alpha * Projection && (displacement_cond) && abs(NewE - previousE) < 10)
    {
      // std::cout<<Bead_1.Total_force.norm()*alpha<<" Displacement of the bead\n";
      break;
    }

    if (std::isnan(E_Vol))
    {
      std::cout << "E vol is nan\n";
    }
    if (std::isnan(E_Sur))
    {
      std::cout << "E sur is nan\n";
    }
    if (std::isnan(E_Ben))
    {
      std::cout << "E ben is nan\n";
    }
    if (std::isnan(E_Bead))
    {
      std::cout << "E bead is nan\n";
    }

    if (std::isnan(NewE))
    {
      std::cout << "New E is nan\n";
      alpha = -1.0;
      break;
    }

    alpha *= rho;

    // if(pulling){

    //   double current_force=Bead_1.Total_force.norm2();
    //   std::cout<<"The force relative dif is "<< abs((current_force-Bead_1.prev_force)/current_force) <<" \n";
    //   std::cout<<"The bead prev force is"<< Bead_1.prev_force<<"\n";
    //   std::cout<<current_force<<"\n";
    //   if(abs((current_force-Bead_1.prev_force)/current_force)<1e-4){
    //   // std::cout<<"\t \t Reseting bead position\n";
    //   // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    //   // // geometry->refreshQuantities();

    //   // X_pos=0.0;
    //   // for(Vertex v : mesh->vertices()){
    //   //   Vertex_pos=geometry->inputVertexPositions[v];
    //   //   if(Vertex_pos.x>X_pos){
    //   //     X_pos=Vertex_pos.x;
    //   //   }
    //   // }
    //   // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
    //   // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){ //|| Bead_1.Pos.x>X_pos+1.4
    //       // Bead_1.strength = Bead_1.strength+0.1;
    //       // std::cout<<"Increasing strength\n";
    //     // }
    //   // Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
    //   // Bead_1.Reset_bead(Vector3({X_pos+1.4,0.0,0.0}));
    //   if(X_pos>40.0){
    //     return -1.0;
    //   }
    // // return alpha;
    // break;
    // }
    // }

    if (alpha < 1e-10)
    {

      // std::cout << "THe timestep got small so the simulation would end \n";
      // std::cout << "THe timestep is " << alpha << " \n";
      // std::cout << "The energy diff is" << abs(NewE - previousE) << "\n";
      // std::cout << "THe relative energy diff  is" << abs((NewE - previousE) / previousE) << "\n";
      // std::cout << "The projection is" << Projection << "\n";
      // std::cout << "The projection is too big " << (Projection > 1.0e10) << " \n";
      if (Projection > 1.0e10)
      {
        // return alpha;
        std::cout << "The gradient got crazy\n";
        std::cout << "The projections is " << Projection << "\n";
        geometry->inputVertexPositions = initial_pos;
        return -1;
      }

      if (!pulling)
      {
        // if(small_TS==true){
        //   alpha=-1.0;

        // }
        if (system_time > 999)
        {
          small_TS = true;
          std::cout << "small timestep\n";
          break;
        }
      }

      // if(Projection>1e10){
      //   std::cout<<"Not moving forward\n";
      //   alpha=0.0;
      // return -1;
      // }
      // continue;
      break;
    }

    geometry->inputVertexPositions = initial_pos + alpha * Force;

    if (pulling && two_bead_pulling)
    {
      // I want to move the beads
      for (size_t i = 0; i < Beads.size(); i++)
      {

        if (Beads[i]->state == "default" || Beads[i]->state == "froze")
        {
          // std::cout<<"TWo here\n";
          Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));
          Beads[i]->Reset_bead(Bead_init[i]);
        }
      }
    }
    // geometry->normalize(Vector3({0.0,0.0,0.0}),false);

    // center = geometry->centerOfMass();

    if (!pulling)
    {
      for (size_t i = 0; i < Beads.size(); i++)
      {
        Beads[i]->Reset_bead(Bead_init[i]);
        Beads[i]->Move_bead(alpha, Vector3({0, 0, 0}));
      }
      // std::cout<<"pulling is false\n";
    }
    // geometry->refreshQuantities();

    A = geometry->totalArea();
    V = geometry->totalVolume();
    D_P = -1 * P0 * (V - V_bar) / (V_bar * V_bar);
    E_Vol = E_Pressure(D_P, V, V_bar);
    E_Sur = E_Surface(KA, A, A_bar);
    E_Ben = E_Bending(H_bar, KB);
    E_Bead = 0;
    for (size_t i = 0; i < Beads.size(); i++)
      E_Bead += Beads[i]->Energy();
    NewE = E_Vol + E_Sur + E_Ben + E_Bead;
    // std::cout<<"THe old energy is "<< previousE <<"\n";
    // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
    // std::cout<<"The projection is :"<<Projection<<"\n";
    // // std::cout<<"THe energy changed to"<<NewE<<"\n";
    // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";
  }

  if (pulling && other_pulling && !two_bead_pulling)
  {
    std::cout << "pulling is true\n";

    //
    double current_force = Beads[0]->Total_force.norm2();
    // double current_force=Bead_1.Total_force.norm2();
    geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
    // geometry->refreshQuantities();
    X_pos = 0.0;
    for (Vertex v : mesh->vertices())
    {
      Vertex_pos = geometry->inputVertexPositions[v];
      if (Vertex_pos.x > X_pos)
      {
        X_pos = Vertex_pos.x;
      }
    }
    // std::cout<<"The system time is"<<system_time<<"\n";
    // if(system_time >2e4){
    std::cout << "The force difference is " << abs((current_force - Bead_1.prev_force) / current_force) << "\n";
    // }
    if (abs((current_force - Beads[0]->prev_force) / current_force) < 1e-5 || alpha < 1e-10)
    {
      // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";

      std::cout << "Some form of steady state was reached\n";
      if (Beads[0]->Pos.x - X_pos < 1.4)
      {
        // Just reset the bead
        std::cout << "\t\t Moving bead\n";
        std::cout << "It is moving " << abs(Beads[0]->Pos.x - (X_pos + 1.4)) << " \n";

        if (abs(Beads[0]->Pos.x - (X_pos + 1.4)) < 1e-6)
        {

          Save_SS = true;
          std::cout << "\t\t Also increasing interaction strength\n";
          Beads[0]->strength = Beads[0]->strength + 0.1;
        }
        Beads[0]->Reset_bead(Vector3({X_pos + 1.4, 0, 0}));
      }
      else
      {
        // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){
        Save_SS = true;
        std::cout << "\t\tIncreasing interaction strength because it is going back\n";
        Beads[0]->strength = Beads[0]->strength + 0.1;
      }

      Beads[0]->prev_E_stationary = E_Bead;

      if (X_pos > 40)
      {
        return -1.0;
      }
      //  Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
      return alpha;
    }

    // Bead_1.Reset_bead(Vector3({X_pos+1.4,0,0}));
  }

  if (pulling && !other_pulling && !two_bead_pulling)
  {
    // std::cout<<"pulling is true\n";

    double current_force = Bead_1.Total_force.norm2();
    // std::cout<<"\t \t Reseting bead position\n";
    geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
    // geometry->refreshQuantities();
    X_pos = 0.0;
    for (Vertex v : mesh->vertices())
    {
      Vertex_pos = geometry->inputVertexPositions[v];
      if (Vertex_pos.x > X_pos)
      {
        X_pos = Vertex_pos.x;
      }
    }
    if (abs((current_force - Bead_1.prev_force) / current_force) < 1e-6 || alpha < 1e-10)
    {
      // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
      std::cout << "It is moving " << abs(Bead_1.Pos.x - (X_pos + 1.4)) << " \n";
      // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){
      std::cout << "Increasing strength\n";
      Save_SS = true;
      std::cout << "The energy difference is " << abs(E_Bead - Bead_1.prev_E_stationary) << " \n";
      if (abs(E_Bead - Bead_1.prev_E_stationary) < 1e-3 && !stop_increasing && Bead_1.strength > 10.0)
      {
        stop_increasing = true;
        std::cout << "\t \t Stopped increasing potential energy\n";
      }
      if (!stop_increasing)
        Bead_1.strength = Bead_1.strength + 0.1;

      Bead_1.prev_E_stationary = E_Bead;

      // }
      Bead_1.Reset_bead(Vector3({X_pos + 1.4, 0, 0}));
      if (X_pos > 40)
      {
        return -1.0;
      }
      //  Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
      return alpha;
    }

    Bead_1.Reset_bead(Vector3({X_pos + 1.4, 0, 0}));
  }

  // std::cout<<"finished while\n";
  if (alpha < 0.0)
  {

    geometry->inputVertexPositions = initial_pos;
  }
  geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);
  // geometry->refreshQuantities();

  return alpha;
}

double Mem3DG::Backtracking_field(VertexData<Vector3> Force, double D_P, double V_bar, double A_bar, double KA, double KB, double H_bar)
{

  // std::cout<<"Backtracking\n";
  bool other_pulling = true;
  // double Bead_force_prev=-12;
  // double Bead_force_new=-15;
  double c1 = 1e-4;
  double rho = 0.5;
  double alpha = 1e-3;
  double positionProjection = 0;
  double X_pos;

  // A=geometry->totalArea();
  // V=geometry->totalVolume();
  // E_Vol = E_Pressure(D_P,V,V_bar);
  // E_Sur = E_Surface(KA,A,A_bar);
  // E_Ben = E_Bending(H_bar,KB);
  // E_Bead = Bead_1.Energy();
  double previousE = E_Vol + E_Sur + E_Ben;
  double NewE;
  VertexData<Vector3> initial_pos(*mesh);

  geometry->normalize(Vector3({0.0, 0.0, 0.0}), false);

  initial_pos = geometry->inputVertexPositions;
  // Vector3 Bead_init = this->Bead_1.Pos;
  // std::cout<<"THe current energy is "<<previousE <<"Is this awful?\n";
  // Zeroth iteration
  double Projection = 0;
  Vector3 center;

  geometry->inputVertexPositions += alpha * Force;
  // geometry->refreshQuantities();

  center = geometry->centerOfMass();
  Vector3 Vertex_pos;

  Total_force = Vector3({0, 0, 0});
  for (Vertex v : mesh->vertices())
  {
    Projection += Force[v.getIndex()].norm2();
    Total_force += Force[v.getIndex()];
  }

  grad_norm = Projection;

  A = geometry->totalArea();
  V = geometry->totalVolume();
  E_Vol = E_Pressure(D_P, V, V_bar);
  E_Sur = E_Surface(KA, A, A_bar);
  E_Ben = E_Bending(H_bar, KB);
  // E_Bead=Bead_1.Energy();
  NewE = E_Vol + E_Sur + E_Ben;

  if (std::isnan(E_Vol))
  {
    std::cout << "E vol is nan\n";
  }
  if (std::isnan(E_Sur))
  {
    std::cout << "E sur is nan\n";
  }
  if (std::isnan(E_Ben))
  {
    std::cout << "E ben is nan\n";
  }

  size_t counter = 0;

  // std::cout<<"starting while\n";
  while (true)
  {
    // if(true){

    if (NewE <= previousE - c1 * alpha * Projection && abs(NewE - previousE) < 10)
    {
      // std::cout<<Bead_1.Total_force.norm()*alpha<<" Displacement of the bead\n";
      break;
    }
    // if(abs(NewE-previousE)>100 && Projection>1e6){
    //   std::cout<<"The energy diff is"<< abs(NewE-previousE)<<"\n";
    //   std::cout<<"THe relative energy diff  is"<<abs((NewE-previousE)/previousE)<<"\n";
    //   std::cout<<"The projection is"<< Projection<<"\n";
    //   std::cout<<"Peak in energy variation, will stop out of safety\n";
    //   alpha=-1;
    //   break;
    // }

    if (std::isnan(E_Vol))
    {
      std::cout << "E vol is nan\n";
    }
    if (std::isnan(E_Sur))
    {
      std::cout << "E sur is nan\n";
    }
    if (std::isnan(E_Ben))
    {
      std::cout << "E ben is nan\n";
    }

    if (std::isnan(NewE))
    {

      alpha = -1.0;
      break;
    }

    alpha *= rho;

    // if(pulling){

    //   double current_force=Bead_1.Total_force.norm2();
    //   std::cout<<"The force relative dif is "<< abs((current_force-Bead_1.prev_force)/current_force) <<" \n";
    //   std::cout<<"The bead prev force is"<< Bead_1.prev_force<<"\n";
    //   std::cout<<current_force<<"\n";
    //   if(abs((current_force-Bead_1.prev_force)/current_force)<1e-4){
    //   // std::cout<<"\t \t Reseting bead position\n";
    //   // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    //   // // geometry->refreshQuantities();

    //   // X_pos=0.0;
    //   // for(Vertex v : mesh->vertices()){
    //   //   Vertex_pos=geometry->inputVertexPositions[v];
    //   //   if(Vertex_pos.x>X_pos){
    //   //     X_pos=Vertex_pos.x;
    //   //   }
    //   // }
    //   // std::cout<<" The difference in distance after a step is "<< abs(Bead_1.Pos.x -(X_pos+1.4)) <<" \n";
    //   // if(abs(Bead_1.Pos.x -(X_pos+1.4))<1e-4 ){ //|| Bead_1.Pos.x>X_pos+1.4
    //       // Bead_1.strength = Bead_1.strength+0.1;
    //       // std::cout<<"Increasing strength\n";
    //     // }
    //   // Bead_1.Reset_bead(Vector3(Bead_1.Pos+Vector3({alpha,0,0})));
    //   // Bead_1.Reset_bead(Vector3({X_pos+1.4,0.0,0.0}));
    //   if(X_pos>40.0){
    //     return -1.0;
    //   }
    // // return alpha;
    // break;
    // }
    // }

    if (alpha < 1e-10)
    {

      // std::cout << "THe timestep got small so the simulation would end \n";
      // std::cout << "THe timestep is " << alpha << " \n";
      // std::cout << "The energy diff is" << abs(NewE - previousE) << "\n";
      // std::cout << "THe relative energy diff  is" << abs((NewE - previousE) / previousE) << "\n";
      // std::cout << "The projection is" << Projection << "\n";

      // if(Projection>1e10){
      //   std::cout<<"Not moving forward\n";
      //   alpha=0.0;
      // return -1;
      // }
      // continue;
      break;
    }

    // geometry->inputVertexPositions = initial_pos+alpha*Force;
    // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
    // // geometry->refreshQuantities();
    // center = geometry->centerOfMass();

    A = geometry->totalArea();
    V = geometry->totalVolume();
    E_Vol = E_Pressure(D_P, V, V_bar);
    E_Sur = E_Surface(KA, A, A_bar);
    E_Ben = E_Bending(H_bar, KB);
    NewE = E_Vol + E_Sur + E_Ben;
    // std::cout<<"THe old energy is "<< previousE <<"\n";
    // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
    // std::cout<<"The projection is :"<<Projection<<"\n";
    // // std::cout<<"THe energy changed to"<<NewE<<"\n";
    // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";
  }

  // std::cout<<"finished while\n";
  if (alpha < 0.0)
  {

    geometry->inputVertexPositions = initial_pos;
  }
  // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  // // geometry->refreshQuantities();

  return alpha;
}

double Mem3DG::Backtracking(VertexData<Vector3> Force, double D_P, double V_bar, double A_bar, double KA, double KB, double H_bar)
{
  double c1 = 1e-4;
  double rho = 0.7;
  double alpha = 1e-3;
  double positionProjection = 0;
  double A = geometry->totalArea();
  double V = geometry->totalVolume();
  double E_Vol = E_Pressure(D_P, V, V_bar);
  double E_Sur = E_Surface(KA, A, A_bar);
  double E_Ben = E_Bending(H_bar, KB);

  // double previousE=E_Vol+E_Sur+E_Ben;
  double previousE = E_Vol + E_Sur + E_Ben;
  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // Zeroth iteration
  double Projection = 0;

  // std::cout<<"THe initial E is "<<previousE<<"\n";
  geometry->inputVertexPositions += alpha * Force;
  // std::cout<< geometry->inputVertexPositions[0]<<"and the other "<< initial_pos[0]<<"\n";

  for (Vertex v : mesh->vertices())
  {
    Projection += Force[v.getIndex()].norm2();
  }

  grad_norm = Projection;

  // geometry->refreshQuantities();

  A = geometry->totalArea();
  V = geometry->totalVolume();
  E_Vol = E_Pressure(D_P, V, V_bar);
  E_Sur = E_Surface(KA, A, A_bar);
  E_Ben = E_Bending(H_bar, KB);

  NewE = E_Vol + E_Sur + E_Ben;
  // NewE=E_Sur+E_Ben;
  // if(std::isnan(E_Vol)){
  //   std::cout<<"E vol is nan\n";
  // }
  if (std::isnan(E_Sur))
  {
    std::cout << "E sur is nan\n";
  }
  if (std::isnan(E_Ben))
  {
    std::cout << "E ben is nan\n";
  }

  size_t counter = 0;
  while (true)
  {
    // if(true){
    // std::cout<<"THe new energy is "<<NewE <<"\n";
    if (NewE <= previousE - c1 * alpha * Projection)
    {
      break;
    }
    // if(std::isnan(E_Vol)){
    // std::cout<<"E vol is nan\n";
    //   }
    if (std::isnan(E_Sur))
    {
      std::cout << "E sur is nan\n";
    }
    if (std::isnan(E_Ben))
    {
      std::cout << "E ben is nan\n";
    }

    alpha *= rho;
    if (alpha < 1e-8)
    {
      // std::cout<<"THe timestep got small\n";
      if (system_time < 2 * Area_evol_steps)
      {
        // std::cout<<"But the area evolution is not complete yet\n";
        break;
      }
      else
      {
        std::cout << "THe simulation will stop because the timestep got smaller than 1e-8 \n";
        alpha = -1.0;
        // continue;
        break;
      }
    }
    // for(Vertex vi : mesh->vertices()){
    //   geometry->inputVertexPositions[vi.getIndex()]= initial_pos[vi.getIndex()]+alpha*Force[vi.getIndex()];
    // }
    geometry->inputVertexPositions = initial_pos + alpha * Force;
    // geometry->refreshQuantities();
    // std::cout<<"THe old energy is "<< previousE <<"\n";
    // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
    // std::cout<<"The projection is :"<<Projection<<"\n";
    // // std::cout<<"THe energy changed to"<<NewE<<"\n";
    // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";

    A = geometry->totalArea();
    V = geometry->totalVolume();
    E_Vol = E_Pressure(D_P, V, V_bar);
    E_Sur = E_Surface(KA, A, A_bar);
    E_Ben = E_Bending(H_bar, KB);
    // NewE=E_Sur+E_Ben;
    NewE = E_Vol + E_Sur + E_Ben;

    if (std::isnan(NewE))
    {
      std::cout << "The energy got Nan\n";

      alpha = -1.0;
      break;
    }
  }

  // for (Edge e : mesh->edges()){
  //   std::cout<< e.remesh;

  // }
  // std::cout<<"\n";
  if (pulling)
  {
    // std::cout<<"THere is pulling right?\n";
    VertexData<Vector3> Horizontal_pull(*mesh, Vector3({1.0, 0.0, 0.0}));

    // std::cout<<"\n";
    geometry->inputVertexPositions += alpha * pulling_force * No_remesh_list_v * Horizontal_pull;
    // geometry->refreshQuantities();
  }
  return alpha;
}

VertexData<Vector3> Mem3DG::Project_force(VertexData<Vector3> Force) const
{

  VertexData<Vector3> Projected_Force(*mesh, Vector3({0, 0, 0}));
  Vector3 Normal_v;
  double projection = 0;
  for (Vertex v : mesh->vertices())
  {
    Normal_v = geometry->vertexNormalAreaWeighted(v);
    projection = dot(Force[v], Normal_v);
    // std::cout<<"The projection is negative?"<< (projection<0) <<" \n";
    Projected_Force[v] = (projection > 0) * projection * Normal_v;
  }
  return Force;
  // return Projected_Force;
}

// THis is backtracking for volume preserving mean curvature flow
double Mem3DG::Backtracking(VertexData<Vector3> Force, double D_P, double V_bar, double KA)
{
  double c1 = 1e-4;
  double rho = 0.7;
  double alpha = 1e-3;
  double positionProjection = 0;
  double A = geometry->totalArea();
  double V = geometry->totalVolume();
  double E_Vol = E_Pressure(D_P, V, V_bar);
  double E_Sur = A * KA;

  // double previousE=E_Vol+E_Sur+E_Ben;
  double previousE = E_Vol + E_Sur;
  double NewE;
  VertexData<Vector3> initial_pos(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // Zeroth iteration
  double Projection = 0;

  // std::cout<<"THe initial E is "<<previousE<<"\n";
  geometry->inputVertexPositions += alpha * Force;
  // std::cout<< geometry->inputVertexPositions[0]<<"and the other "<< initial_pos[0]<<"\n";

  for (Vertex v : mesh->vertices())
  {
    Projection += Force[v.getIndex()].norm2();
  }

  grad_norm = Projection;

  // geometry->refreshQuantities();

  A = geometry->totalArea();
  V = geometry->totalVolume();
  E_Vol = E_Pressure(D_P, V, V_bar);
  E_Sur = A * KA;
  NewE = E_Vol + E_Sur;
  // NewE=E_Sur+E_Ben;
  // if(std::isnan(E_Vol)){
  //   std::cout<<"E vol is nan\n";
  // }
  if (std::isnan(E_Sur))
  {
    std::cout << "E sur is nan\n";
  }

  size_t counter = 0;
  while (true)
  {
    // if(true){
    // std::cout<<"THe new energy is "<<NewE <<"\n";
    if (NewE <= previousE - c1 * alpha * Projection)
    {
      break;
    }
    // if(std::isnan(E_Vol)){
    // std::cout<<"E vol is nan\n";
    //   }
    if (std::isnan(E_Sur))
    {
      std::cout << "E sur is nan\n";
    }

    if (std::isnan(NewE))
    {
      std::cout << "The energy got Nan\n";

      alpha = -1.0;
      break;
    }

    alpha *= rho;
    if (alpha < 1e-8)
    {
      // std::cout<<"THe timestep got small\n";
      if (system_time < 2 * Area_evol_steps)
      {
        // std::cout<<"But the area evolution is not complete yet\n";
        break;
      }
      else
      {
        std::cout << "THe simulation will stop because the timestep got smaller than 1e-8 \n";
        alpha = -1.0;
        // continue;
        break;
      }
    }
    // for(Vertex vi : mesh->vertices()){
    //   geometry->inputVertexPositions[vi.getIndex()]= initial_pos[vi.getIndex()]+alpha*Force[vi.getIndex()];
    // }
    geometry->inputVertexPositions = initial_pos + alpha * Force;
    // geometry->refreshQuantities();
    // std::cout<<"THe old energy is "<< previousE <<"\n";
    // std::cout<<"Alpha is "<< alpha<<"and the new energy is"<< NewE << "\n";
    // std::cout<<"The projection is :"<<Projection<<"\n";
    // // std::cout<<"THe energy changed to"<<NewE<<"\n";
    // std::cout<< "Volume E"<<E_Vol <<"Surface E" << E_Sur <<"\n";

    A = geometry->totalArea();
    V = geometry->totalVolume();
    E_Vol = E_Pressure(D_P, V, V_bar);
    E_Sur = A * KA;
    NewE = E_Vol + E_Sur;
  }

  // for (Edge e : mesh->edges()){
  //   std::cout<< e.remesh;

  // }
  // std::cout<<"\n";
  if (pulling)
  {
    // std::cout<<"THere is pulling right?\n";
    VertexData<Vector3> Horizontal_pull(*mesh, Vector3({1.0, 0.0, 0.0}));

    // std::cout<<"\n";
    geometry->inputVertexPositions += alpha * pulling_force * No_remesh_list_v * Horizontal_pull;
    // geometry->refreshQuantities();
  }
  return alpha;
}

double Mem3DG::integrate_Newton_Normal_Sherman(std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data, std::vector<std::string> Constraints, std::vector<std::string> Data_filenames)
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

  if (Newton_iter == 0)
  {
    // First newton iter
    mu = 1e-2;
    // This mu is the value for the area and volume constraint
  }
  // Here i will update the values for the area and volume constraint c:
  for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
  {
    if (Sim_handler->Energies[i] == "Volume_constraint")
    {
      Sim_handler->Energy_constants[i][0] = mu;
    }
    if (Sim_handler->Energies[i] == "Area_constraint")
    {
      Sim_handler->Energy_constants[i][0] = mu;
    }
  }
  // SO now the volume and area constraints are back back in the game.

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
  {
    // We need to add the value of the constraint
    if (Constraints[Ci] == "Volume")
    {
      RHS((N_vert) + Ci) = -1 * (geometry->totalVolume() - Sim_handler->Trgt_vol);
      // RHS((N_vert + N_beads) + Ci) = -1 * (geometry->totalVolume() - Sim_handler->Trgt_vol);

      std::cout << "The total volume is " << geometry->totalVolume() << " and the target is " << Sim_handler->Trgt_vol << "\n";
    }
    if (Constraints[Ci] == "Area")
    {
      RHS((N_vert) + Ci) = -1 * (A - Sim_handler->Trgt_area);

      std::cout << "The total area  is " << A << " and the target is " << Sim_handler->Trgt_area << "\n";
    }
    if (Constraints[Ci] == "CMx")
      RHS((N_vert) + Ci) = 0.0; // Considering the beads should have +N_beads
    if (Constraints[Ci] == "CMy")
      RHS((N_vert) + Ci) = 0.0;
    if (Constraints[Ci] == "CMz")
      RHS((N_vert) + Ci) = 0.0;
    if (Constraints[Ci] == "Rx")
      RHS((N_vert) + Ci) = 0.0;
    if (Constraints[Ci] == "Ry")
      RHS((N_vert) + Ci) = 0.0;
    if (Constraints[Ci] == "Rz")
      RHS((N_vert) + Ci) = 0.0;
  }
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

  //

  double Projection = result.transpose() * LHS * RHS;
  double Current_grad_norm = 0.0;

  Sim_handler->Calculate_Lag_norm_Normal(&Current_grad_norm);

  Sim_handler->Current_grad = Force_result;

  // Current grad norm is

  double backtrackstep;
  // if (false)
  // {
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
      std::cout << "Projection is " << Projection << " and the grad norm is " << Current_grad_norm << "\n";
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
    A = 0.0;
    geometry->requireFaceAreas();
    for (Face f : mesh->faces())
    {
      A += geometry->faceArea(f);
    }
    Sim_data << time + backtrackstep << " " << discreteTs << " " << geometry->totalVolume() << " " << A << " ";
    for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
    {

      Sim_data << Sim_handler->Energy_values[i] << " ";
    }

    Sim_data << TotE << " ";

    // Sim_data << Current_grad_norm <<" " << backtrackstep <<" \n";
    for (size_t i = 0; i < Sim_handler->Gradient_norms.size(); i++)
    {
      Sim_data << Sim_handler->Gradient_norms[i] << " ";
    }
    // Sim_data << Grad_tot_norm << " ";
    Sim_data << backtrackstep << " \n";
  }

  if (Bead_data_filenames.size() != 0 && Save_output_data)
  {

    std::ofstream Bead_data;
    for (size_t i = 0; i < Beads.size(); i++)
    {
      Bead_data = std::ofstream(Bead_data_filenames[i], std::ios_base::app);
      Bead_data << discreteTs << " " << Beads[i]->Pos.x << " " << Beads[i]->Pos.y << " " << Beads[i]->Pos.z << " " << Beads[i]->Total_force.x << " " << Beads[i]->Total_force.y << " " << Beads[i]->Total_force.z << " \n";
      Bead_data.close();
    }
  }

  return backtrackstep;
  // Ok so we have everything i guess?
}

double Mem3DG::integrate_implicit(std::vector<std::string> Energies, std::vector<std::vector<double>> Energy_constants, std::ofstream &Sim_data, double time, std::vector<std::string> Bead_data_filenames, bool Save_output_data)
{

  // Okok what now

  return 0.0;
}

/*
 * Performs integration with the bead
 *
 * Input: The timestep <h>.
 * Returns:
 */
double Mem3DG::integrate(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd, std::ofstream &Sim_data, double time, bool bead, std::vector<std::string> Bead_data_filenames, bool Save_output_data, bool pulling)
{
  //, Beads Bead_1
  std::ofstream Bead_data;
  // std::cout<<"1_2\n";
  Bead *Active_bead;
  if (bead)
  {
    // std::cout<<"This is being called\n";
    // std::cout<<"There are "<< Beads.size() <<" beads\n";
    for (size_t i = 0; i < Beads.size(); i++)
    {

      Bead_data = std::ofstream(Bead_data_filenames[i], std::ios_base::app);
      Bead_data << Beads[i]->Pos.x << " " << Beads[i]->Pos.y << " " << Beads[i]->Pos.z << " " << Beads[i]->Total_force.x << " " << Beads[i]->Total_force.y << " " << Beads[i]->Total_force.z << " \n";
      // std::cout<<Bead_1.Pos.x << " "<< Bead_1.Pos.y << " "<< Bead_1.Pos.z<<" \n";
      // std::cout<<"The total force is "<<Bead_1.Total_force <<"\n";
      Bead_data.close();
    }
  }

  // std::cout<<"2_2\n";
  // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);

  VertexData<Vector3> Force(*mesh);

  Force = buildFlowOperator(h, V_bar, nu, c0, P0, KA, KB, Kd); //+Bead_1.Gradient();

  // std::cout<<"3_2\n";
  // std::cout<<"\t \tThe size of Beads is"<< Beads.size()<<"\n";
  VertexData<Vector3> Bead_force;
  for (size_t i = 0; i < Beads.size(); i++)
  {
    // std::cout<<"Iterating over bead number "<<i+1<<"\n";
    Bead_force = Beads[i]->Gradient();
    Beads[i]->Bead_interactions();
    Force += Bead_force;
  }
  // VertexData<Vector3> Bead_force = Bead_1.Gradient();

  // VertexData<Vector3> Bead_force=Project_force(Bead_1.Gradient());

  // std::cout<<"4_2\n";
  // Force+=Bead_force;

  // This force is the grdient basically
  double alpha = h;

  // I have the forces for almost everything i just need the bead.

  size_t vindex;
  size_t Nvert = mesh->nVertices();

  // I want to print the Volume, the area, VOl_E, Area_E, Bending_E
  V = geometry->totalVolume();
  A = geometry->totalArea();
  double A_bar = 4 * PI * pow(3 * V_bar / (4 * PI * nu), 2.0 / 3.0);
  double H_bar = sqrt(4 * PI / A_bar) * c0 / 2.0; // Coment this with another comment
  double D_P = -1 * P0 * (V - V_bar) / (V_bar * V_bar);
  double lambda = KA * (A - A_bar) / (A_bar * A_bar);

  E_Vol = E_Pressure(D_P, V, V_bar);
  E_Sur = E_Surface(KA, A, A_bar);
  E_Ben = E_Bending(H_bar, KB);
  E_Bead = 0.0;

  for (size_t i = 0; i < Beads.size(); i++)
  {

    E_Bead += Beads[i]->Energy();
  }

  double backtrackstep;

  // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";

  // if(system_time< 50){
  //   backtrackstep=h;
  // }
  // else{
  // std::cout<<"Backtracking is being called\n";
  Total_force = Vector3({0, 0, 0});
  backtrackstep = Backtracking(Force, P0, V_bar, A_bar, KA, KB, H_bar, bead, pulling);
  // Total_force*=backtrackstep;
  // I need to find the stop increasing bt
  // }
  // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";

  // geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;

  // Bead_1.Move_bead(backtrackstep,center);

  if ((Save_output_data) || backtrackstep < 0)
  {
    Sim_data << V_bar << " " << A_bar << " " << time << " " << discreteTs << " " << V << " " << A << " " << E_Vol << " " << E_Sur << " " << E_Ben << " " << E_Bead << " " << grad_norm << " " << backtrackstep << " \n";
  }
  system_time += 1;

  return backtrackstep;
}

double Mem3DG::integrate_field(double h, double V_bar, double nu, double P0, double KA, double KB, double slope, double x0, std::ofstream &Sim_data, double time, bool Save)
{
  //, Beads Bead_1

  // std::cout<<"2_2\n";
  // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
  VertexData<Vector3> Force(*mesh);
  VertexData<Vector3> Field(*mesh);
  double c0 = 0.0;
  double Kd = 0.0;
  Force = buildFlowOperator(V_bar, P0, KA, KB, h); //+Bead_1.Gradient();

  // This force is the grdient basically
  double alpha = h;

  // I have the forces for almost everything i just need the bead.

  size_t vindex;
  size_t Nvert = mesh->nVertices();

  // I want to print the Volume, the area, VOl_E, Area_E, Bending_E
  V = geometry->totalVolume();
  A = geometry->totalArea();
  double A_bar = 4 * PI * pow(3 * V_bar / (4 * PI * nu), 2.0 / 3.0);
  double H_bar = sqrt(4 * PI / A_bar) * c0 / 2.0; // Coment this with another comment
  double D_P = -1 * P0 * (V - V_bar) / (V_bar * V_bar);
  double lambda = KA * (A - A_bar) / (A_bar * A_bar);

  E_Vol = E_Pressure(D_P, V, V_bar);
  E_Sur = E_Surface(KA, A, A_bar);
  E_Ben = E_Bending(H_bar, KB);
  E_Bead = Bead_1.Energy();

  // double backtrackstep;

  // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";

  // if(system_time< 50){
  //   backtrackstep=h;
  // }
  // else{
  // std::cout<<"Backtracking is being called\n";
  Total_force = Vector3({0, 0, 0});
  // backtrackstep=Backtracking(Force,D_P,V_bar,A_bar,KA,KB,H_bar,bead,pulling);
  // Total_force*=backtrackstep;
  // I need to find the stop increasing bt
  // }
  // std::cout<<" The position of vertex 120 is "<<geometry->inputVertexPositions[120].x<<" "<<geometry->inputVertexPositions[120].y<< " "<<geometry->inputVertexPositions[120].z<<" \n";
  // double backtrackstep = Backtracking_field(Force,D_P,V_bar,A_bar,KA,KB,H_bar);
  double backtrackstep = 1e-4; // lowest 5e-7
  // I need to calculate the force field

  Field = Linear_force_field(x0, slope);
  // I want to make it volume preserving

  Force += Field;
  geometry->inputVertexPositions = geometry->inputVertexPositions + backtrackstep * Force;

  // I want to find the leftmost particle

  Vector3 Leftmost({0, 0, 0});
  Vector3 Rightmost({0, 0, 0});

  int rightmost;
  int leftmost;
  for (Vertex v : mesh->vertices())
  {
    if (geometry->inputVertexPositions[v].x < Leftmost.x)
    {
      Leftmost = geometry->inputVertexPositions[v];
      leftmost = v.getIndex();
    }
    if (geometry->inputVertexPositions[v].x > Rightmost.x)
    {
      Rightmost = geometry->inputVertexPositions[v];
      rightmost = v.getIndex();
    }
  }
  Leftmost -= Vector3({-5, 0, 0});
  VertexData<Vector3> Displacement(*mesh, Leftmost);

  geometry->inputVertexPositions -= Displacement;

  // std::cout<<"The force at the leftmost is"<< Field[leftmost] <<" with magnitude "<< Field[leftmost].norm()<<"\n";
  // std::cout<<"The force at the leftmost is"<< Field[rightmost] <<" with magnitude "<< Field[rightmost].norm()<<"\n";

  // std::cout<<" The difference in force in theory " << slope*(geometry->inputVertexPositions[rightmost].x - geometry->inputVertexPositions[leftmost].x)<< " and the obtained is " << Field[rightmost].norm()-Field[leftmost].norm()<<"\n";

  // geometry->normalize(Vector3({0.0,0.0,0.0}),false);
  // geometry->refreshQuantities();

  // Bead_1.Move_bead(backtrackstep,center);

  if (Save)
  {
    Sim_data << V_bar << " " << A_bar << " " << time << " " << V << " " << A << " " << E_Vol << " " << E_Sur << " " << E_Ben << " " << E_Bead << " " << grad_norm << " " << backtrackstep << " \n";
  }
  system_time += 1;

  return backtrackstep;
}

// THis is the main integrate
double Mem3DG::integrate(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd, std::ofstream &Sim_data, double time, bool Save)
{

  VertexData<Vector3> Force(*mesh);
  Force = buildFlowOperator(h, V_bar, nu, c0, P0, KA, KB, Kd);

  // This force is the gradient basically
  double alpha = 1e-3;

  // I have the forces for almost everything i just need the bead.

  size_t vindex;
  size_t Nvert = mesh->nVertices();

  // I want to print the Volume, the area, VOl_E, Area_E, Bending_E
  double V = geometry->totalVolume();
  double A = geometry->totalArea();
  double A_bar = 4 * PI * pow(3 * V_bar / (4 * PI * nu), 2.0 / 3.0);
  double H_bar = sqrt(4 * PI / A_bar) * c0 / 2.0; // Coment this with another comment
  double D_P = -1 * P0 * (V - V_bar) / V_bar / V_bar;
  double lambda = KA * (A - A_bar) / A_bar;

  double E_Vol = E_Pressure(D_P, V, V_bar);
  double E_Sur = E_Surface(KA, A, A_bar);
  double E_Ben = E_Bending(H_bar, KB);
  double backtrackstep;

  // if(system_time< 50){
  //   backtrackstep=h;
  // }
  // else{
  backtrackstep = Backtracking(Force, D_P, V_bar, A_bar, KA, KB, H_bar);
  // }
  if (Save || backtrackstep < 0)
  {
    Sim_data << V_bar << " " << A_bar << " " << time << " " << V << " " << A << " " << E_Vol << " " << E_Sur << " " << E_Ben << " " << grad_norm << " " << backtrackstep << " \n";
  }
  system_time += 1;

  // for (Vertex v : mesh->vertices()) {
  //     vindex=v.getIndex();

  //     // Update= { Total_force[vindex],Total_force[vindex+Nvert],Total_force[vindex+2*Nvert]  };
  //     // Update=h*Force[vindex];
  //     Update=backtrackstep*Force[vindex];
  //     geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
  // }
  // Force=MeshData<Vertex,Vector3>::fromVector(Delta_x);
  // geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;

  return backtrackstep;
}

// THis integrate will do surface tension + volume preservation
double Mem3DG::integrate(double h, double V_bar, double P0, double KA, std::ofstream &Sim_data, double time, bool Save)
{

  VertexData<Vector3> Force(*mesh);
  Force = buildFlowOperator(h, V_bar, P0, KA);

  // This force is the gradient basically
  double alpha = 1e-2;

  // I have the forces for almost everything i just need the bead.

  size_t vindex;
  size_t Nvert = mesh->nVertices();

  // I want to print the Volume, the area, VOl_E, Area_E, Bending_E
  double V = geometry->totalVolume();
  double A = geometry->totalArea();
  double D_P = -1 * P0 * (V - V_bar) / V_bar / V_bar;

  double E_Vol = E_Pressure(D_P, V, V_bar);
  double E_Sur = KA * A;
  double backtrackstep;

  backtrackstep = Backtracking(Force, D_P, V_bar, KA);
  // }
  if (Save || backtrackstep < 0)
  {
    Sim_data << V_bar << " " << time << " " << V << " " << A << " " << E_Vol << " " << E_Sur << " " << grad_norm << " " << backtrackstep << " \n";
  }
  system_time += 1;

  // for (Vertex v : mesh->vertices()) {
  //     vindex=v.getIndex();

  //     // Update= { Total_force[vindex],Total_force[vindex+Nvert],Total_force[vindex+2*Nvert]  };
  //     // Update=h*Force[vindex];
  //     Update=backtrackstep*Force[vindex];
  //     geometry->inputVertexPositions[v] =geometry->inputVertexPositions[v]+ Update ; // placeholder
  // }
  // Force=MeshData<Vertex,Vector3>::fromVector(Delta_x);
  // geometry->inputVertexPositions=geometry->inputVertexPositions+backtrackstep*Force;

  return backtrackstep;
}

void Mem3DG::Grad_Vol_dx(std::ofstream &Gradient_file, double P0, double V_bar, size_t index) const
{
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // VertexData<Vector3> Gradients(*mesh,0.0);
  double V = geometry->totalVolume();
  double E_vol = 0.0;
  double E_vol_back = 0.0;
  double E_vol_front = 0.0;
  double dr;
  Vector3 grad{0.0, 0.0, 0.0};
  double D_P = -P0 * (V - V_bar) / (V_bar * V_bar);
  E_vol = E_Pressure(D_P, V, V_bar);

  VertexData<Vector3> Calc_grad = D_P * OsmoticPressure();

  Gradient_file << Calc_grad[index].x << " " << Calc_grad[index].y << " " << Calc_grad[index].z << " \n";

  Vector<Vector3> Gradients(10);
  dr = 10.0;
  for (size_t exponent = 0; exponent < 20; exponent++)
  {
    dr = dr / 10.0;
    // dr=pow(10,-1*exponent);

    // std::cout<<dr<<" ";

    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_front = E_Pressure(D_P, V, V_bar);
    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_back = E_Pressure(D_P, V, V_bar);
    grad.x = (E_vol_front - E_vol_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_back = E_Pressure(D_P, V, V_bar);
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_front = E_Pressure(D_P, V, V_bar);
    grad.y = (E_vol_front - E_vol_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_back = E_Pressure(D_P, V, V_bar);
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_front = E_Pressure(D_P, V, V_bar);
    grad.z = (E_vol_front - E_vol_back) / (2 * dr);

    Gradients[exponent] = grad;
    Gradient_file << Gradients[exponent].x << " " << Gradients[exponent].y << " " << Gradients[exponent].z << " \n";
  }

  // Gradient_file<<"\n";
  // This function should just write
  // std::cout<<"\n";

  return;
}

VertexData<Vector3> Mem3DG::Grad_Vol(std::ofstream &Gradient_file, double P0, double V_bar, bool Save) const
{
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  VertexData<Vector3> Finite_grad(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // VertexData<Vector3> Gradients(*mesh,0.0);
  double V = geometry->totalVolume();
  double E_vol = 0.0;
  double E_vol_back = 0.0;
  double E_vol_front = 0.0;
  double dr;
  double D_P = -P0 * (V - V_bar) / V_bar / V_bar;
  double total_grad_finite = 0;
  double total_grad_theory = 0;

  Vector3 grad{0.0, 0.0, 0.0};
  Vector3 grad_theory;
  Vector3 difference;

  E_vol = E_Pressure(D_P, V, V_bar);

  VertexData<Vector3> Calc_grad = D_P * OsmoticPressure();

  dr = 1e-6;

  size_t N_vert = mesh->nVertices();
  // for(size_t index=0; index<N_vert; index++){
  size_t index;
  for (Vertex v : mesh->vertices())
  {
    index = v.getIndex();
    grad_theory = Calc_grad[v];
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_front = E_Pressure(D_P, V, V_bar);
    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_back = E_Pressure(D_P, V, V_bar);
    grad.x = (E_vol_front - E_vol_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_back = E_Pressure(D_P, V, V_bar);
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_front = E_Pressure(D_P, V, V_bar);
    grad.y = (E_vol_front - E_vol_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_back = E_Pressure(D_P, V, V_bar);
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    V = geometry->totalVolume();
    D_P = -P0 * (V - V_bar) / V_bar / V_bar;
    E_vol_front = E_Pressure(D_P, V, V_bar);
    grad.z = (E_vol_front - E_vol_back) / (2 * dr);
    if (Save)
    {
      difference = grad + grad_theory;
      Gradient_file << difference.x << " " << difference.y << " " << difference.z << " " << difference.norm() / grad.norm() << " " << grad.norm() / grad_theory.norm() << " \n";
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
    Gradient_file << sqrt(total_grad_theory) << " " << sqrt(total_grad_finite) << "\n";
  }

  return Finite_grad;
}

VertexData<Vector3> Mem3DG::Grad_Area(std::ofstream &Gradient_file, double A_bar, double KA, bool Save) const
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
  double lambda = KA * (A - A_bar) / A_bar;
  double total_grad_finite = 0;
  double total_grad_theory = 0;

  Vector3 grad{0.0, 0.0, 0.0};
  Vector3 grad_theory;
  Vector3 difference;

  E_area = E_Surface(KA, A, A_bar);

  VertexData<Vector3> Calc_grad = lambda * SurfaceTension();

  dr = 1e-6;

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
    lambda = KA * (A - A_bar) / A_bar;
    E_area_front = E_Surface(KA, A, A_bar);
    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    lambda = KA * (A - A_bar) / A_bar;
    E_area_back = E_Surface(KA, A, A_bar);
    grad.x = (E_area_front - E_area_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    lambda = KA * (A - A_bar) / A_bar;
    E_area_back = E_Surface(KA, A, A_bar);
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    lambda = KA * (A - A_bar) / A_bar;
    E_area_front = E_Surface(KA, A, A_bar);
    grad.y = (E_area_front - E_area_back) / (2 * dr);

    geometry->inputVertexPositions[v] = initial_pos[v] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    lambda = KA * (A - A_bar) / A_bar;
    E_area_back = E_Surface(KA, A, A_bar);
    geometry->inputVertexPositions[v] = initial_pos[v] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    A = geometry->totalArea();
    lambda = KA * (A - A_bar) / A_bar;
    E_area_front = E_Surface(KA, A, A_bar);
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
    Gradient_file << sqrt(total_grad_theory) << " " << sqrt(total_grad_finite) << "\n";
  }

  return Finite_grad;
}

VertexData<Vector3> Mem3DG::Grad_Bending(std::ofstream &Gradient_file, double H_bar, double KB, bool Save)
{
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  VertexData<Vector3> Finite_grad(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // VertexData<Vector3> Gradients(*mesh,0.0);
  double E_bending = 0.0;
  double E_bending_back = 0.0;
  double E_bending_front = 0.0;
  double dr;

  Vector3 grad{0.0, 0.0, 0.0};
  Vector3 grad_theory;
  Vector3 difference;

  E_bending = E_Bending(H_bar, KB);

  double norm_whole_finite = 0;
  double norm_whole_theory = 0;

  VertexData<Vector3> Calc_grad = KB * Bending(H_bar);

  dr = 1e-6;

  size_t N_vert = mesh->nVertices();

  for (size_t index = 0; index < N_vert; index++)
  {
    grad_theory = Calc_grad[index];
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();

    E_bending_front = E_Bending(H_bar, KB);
    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();

    E_bending_back = E_Bending(H_bar, KB);
    grad.x = (E_bending_front - E_bending_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();

    E_bending_back = E_Bending(H_bar, KB);
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();

    E_bending_front = E_Bending(H_bar, KB);
    grad.y = (E_bending_front - E_bending_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();

    E_bending_back = E_Bending(H_bar, KB);
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();

    E_bending_front = E_Bending(H_bar, KB);
    grad.z = (E_bending_front - E_bending_back) / (2 * dr);
    if (Save)
    {
      difference = grad + grad_theory;
      Gradient_file << difference.x << " " << difference.y << " " << difference.z << " " << difference.norm() / grad.norm() << " " << grad.norm() / grad_theory.norm() << " \n";
      // difference= grad_theory;
      // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" "<<grad_theory.norm()<<" \n" ;
      norm_whole_finite += grad.norm2();
      norm_whole_theory += grad_theory.norm2();
    }
    Finite_grad[index] = -1 * grad;
    geometry->inputVertexPositions[index] = initial_pos[index];
    // geometry->refreshQuantities();
  }

  norm_whole_finite = sqrt(norm_whole_finite);
  norm_whole_theory = sqrt(norm_whole_theory);
  if (Save)
  {
    Gradient_file << norm_whole_finite << " " << norm_whole_finite << " \n";
  }

  return Finite_grad;
}

void Mem3DG::Grad_Bending_2(std::ofstream &Gradient_file, double H_bar, double KB)
{
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // VertexData<Vector3> Gradients(*mesh,0.0);
  double E_bending = 0.0;
  double E_bending_back = 0.0;
  double E_bending_front = 0.0;
  double dr;

  Vector3 grad{0.0, 0.0, 0.0};
  Vector3 grad_theory;
  Vector3 difference;

  E_bending = E_Bending(H_bar, KB);

  double norm_diff = 0;

  VertexData<Vector3> Calc_grad = KB * Bending(H_bar);

  dr = 1e-6;

  size_t N_vert = mesh->nVertices();

  for (size_t index = 0; index < N_vert; index++)
  {
    grad_theory = Calc_grad[index];
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();

    E_bending_front = E_Bending(H_bar, KB);
    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();

    E_bending_back = E_Bending(H_bar, KB);
    grad.x = (E_bending_front - E_bending_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();

    E_bending_back = E_Bending(H_bar, KB);
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();

    E_bending_front = E_Bending(H_bar, KB);
    grad.y = (E_bending_front - E_bending_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();

    E_bending_back = E_Bending(H_bar, KB);
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();

    E_bending_front = E_Bending(H_bar, KB);
    grad.z = (E_bending_front - E_bending_back) / (2 * dr);

    difference = grad + grad_theory;

    // Gradient_file<< difference.x <<" "<<difference.y<<" "<< difference.z<<" \n" ;
    // difference= grad_theory;
    Gradient_file << difference.x << " " << difference.y << " " << difference.z << " " << difference.norm() / grad_theory.norm() << " \n";
    // norm_whole_finite+=grad.norm2();
    // norm_whole_theory+=grad_theory.norm2();
  }

  // norm_whole_finite=sqrt(norm_whole_finite);
  // norm_whole_theory=sqrt(norm_whole_theory);

  // Gradient_file<< norm_whole_finite<<" "<< norm_whole_finite <<" \n";

  return;
}

void Mem3DG::Bending_test(std::ofstream &Analysis_file, double H0, double KB)
{

  // THe idea here would be to do the same as with the other function but store things

  size_t neigh_index;
  size_t N_vert = mesh->nVertices();

  Vector3 Hij;
  Vector3 Kij;
  Vector3 Sij_1;
  Vector3 Sij_2;
  Vector3 Sij_22;

  Vector3 F1 = {0, 0, 0};
  Vector3 F2 = {0, 0, 0};
  Vector3 F3 = {0, 0, 0};
  Vector3 F4 = {0, 0, 0};

  Vector3 Position_1;
  Vector3 Position_2;

  VertexData<Vector3> Force(*mesh);

  VertexData<double> Scalar_MC(*mesh, 0.0);
  double factor;
  size_t index1;
  for (Vertex v1 : mesh->vertices())
  {
    index1 = v1.getIndex();
    Scalar_MC[index1] = geometry->scalarMeanCurvature(v1) / geometry->barycentricDualArea(v1);
  }

  size_t index;
  for (Vertex v : mesh->vertices())
  {
    F1 = {0, 0, 0};
    F2 = {0, 0, 0};
    F3 = {0, 0, 0};
    F4 = {0, 0, 0};

    index = v.getIndex();
    Position_1 = geometry->inputVertexPositions[v];
    for (Halfedge he : v.outgoingHalfedges())
    {

      neigh_index = he.tipVertex().getIndex();
      Position_2 = geometry->inputVertexPositions[neigh_index];

      Kij = computeHalfedgeGaussianCurvatureVector(he);
      // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))-1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
      factor = -1 * (Scalar_MC[index] - H0) - 1 * (Scalar_MC[neigh_index] - H0);
      F1 = F1 + factor * Kij;

      Hij = 2 * computeHalfedgeMeanCurvatureVector(he);
      // factor=(1/3.0)*(Scalar_MC[index]*Scalar_MC[index] -(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0)*(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0))+(2.0/3.0)*(Scalar_MC[neigh_index]*Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0)*(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
      factor = (1 / 3.0) * (Scalar_MC[index] * Scalar_MC[index] - H0 * H0) + (2.0 / 3.0) * (Scalar_MC[neigh_index] * Scalar_MC[neigh_index] - H0 * H0);
      F2 = F2 + factor * Hij;

      Sij_1 = geometry->edgeLength(he.edge()) * dihedralAngleGradient(he, he.vertex());
      // factor=-1*(Scalar_MC[index]-(system_time<50? H_Vector_0[index]+dH_Vector[index]*system_time: H0));
      factor = -1 * (Scalar_MC[index] - H0);
      F3 = F3 + factor * Sij_1;

      Sij_2 = (geometry->edgeLength(he.twin().edge()) * dihedralAngleGradient(he.twin(), he.vertex()) + geometry->edgeLength(he.next().edge()) * dihedralAngleGradient(he.next(), he.vertex()) + geometry->edgeLength(he.twin().next().next().edge()) * dihedralAngleGradient(he.twin().next().next(), he.vertex()));

      // Sij_2=-1*( geometry->cotan(he.next().next())*geometry->faceNormal(he.face()) + geometry->cotan(he.twin())*geometry->faceNormal(he.twin().face()));

      // std::cout<<"Norm new "<< Sij_2.norm()<<"Norm old "<< Sij_22.norm()<<"\n";
      // std::cout<<"Norm ratio new/old = "<<Sij_2.norm()/Sij_22.norm()<<" \n";
      // std::cout<<"Relative orientation "<<dot(Sij_2/Sij_2.norm(),Sij_22/Sij_22.norm())<<" \n\n";

      // factor= -1*(Scalar_MC[neigh_index]-(system_time<50? H_Vector_0[neigh_index]+dH_Vector[neigh_index]*system_time: H0));
      factor = -1 * (Scalar_MC[neigh_index] - H0);

      F4 = F4 + factor * Sij_2;
    }
    Analysis_file << F1.x << " " << F1.y << " " << F1.z << " " << F2.x << " " << F2.y << " " << F2.z << " " << F3.x << " " << F3.y << " " << F3.z << " " << F4.x << " " << F4.y << " " << F4.z << " \n";

    // Force[index]=F1+F2+F3+F4;
  }

  return;
}

double Mem3DG::integrate_finite(double h, double V_bar, double nu, double c0, double P0, double KA, double KB, double Kd, std::ofstream &Sim_data, double time, std::ofstream &Gradient_file_vol, std::ofstream &Gradient_file_area, std::ofstream &Gradient_file_bending, bool Save)
{

  // Vector<double> Total_force=buildFlowOperator(h,V_bar,nu,c0,P0,KA,KB,Kd);
  VertexData<Vector3> Force(*mesh);
  Vector3 Update;
  size_t vindex;
  size_t Nvert = mesh->nVertices();

  // I want to print the Volume, the area, VOl_E, Area_E, Bending_E
  double V = geometry->totalVolume();
  double A = geometry->totalArea();
  double A_bar = 4 * PI * pow(3 * V_bar / (4 * PI * nu), 2.0 / 3.0);
  double H_bar = sqrt(4 * PI / A_bar) * c0 / 2.0; // Coment this with another comment
  double D_P = -1 * P0 * (V - V_bar) / V_bar / V_bar;
  double lambda = KA * (A - A_bar) / A_bar;
  // std::cout<<"The forces are calculated here\n";
  Force = Grad_Vol(Gradient_file_vol, P0, V_bar, Save) + Grad_Area(Gradient_file_area, A_bar, KA, Save) + Grad_Bending(Gradient_file_bending, H_bar, KB, Save);
  // lambda=KA;
  // KB=0;

  double E_Vol = E_Pressure(D_P, V, V_bar);
  double E_Sur = E_Surface(KA, A, A_bar);
  double E_Ben = E_Bending(H_bar, KB);
  double backtrackstep;

  // std::cout<<"Now is the backtracking\n";

  backtrackstep = Backtracking(Force, D_P, V_bar, A_bar, KA, KB, H_bar);
  Sim_data << V_bar << " " << A_bar << " " << time << " " << V << " " << A << " " << E_Vol << " " << E_Sur << " " << E_Ben << " " << grad_norm << " " << backtrackstep << " " << "\n";

  system_time += 1;

  for (Vertex v : mesh->vertices())
  {
    vindex = v.getIndex();
    Update = backtrackstep * Force[vindex];
    geometry->inputVertexPositions[v] = geometry->inputVertexPositions[v] + Update; // placeholder
  }

  return backtrackstep;
}

void Mem3DG::Grad_Bead_dx(std::ofstream &Gradient_file, bool Save)
{
  // I want to calculate the gradient of the volume
  VertexData<Vector3> initial_pos(*mesh);
  initial_pos = geometry->inputVertexPositions;
  // VertexData<Vector3> Gradients(*mesh,0.0);
  // double V=geometry->totalVolume();
  double E_bead = 0.0;
  double E_bead_back = 0.0;
  double E_bead_front = 0.0;
  double dr;
  Vector3 grad{0.0, 0.0, 0.0};
  E_bead = Bead_1.Energy();
  int index = 7333;

  VertexData<Vector3> Calc_grad = Bead_1.Gradient();

  Gradient_file << Calc_grad[index].x << " " << Calc_grad[index].y << " " << Calc_grad[index].z << " \n";

  Vector<Vector3> Gradients(20);
  dr = 10.0;
  for (size_t exponent = 0; exponent < 20; exponent++)
  {
    dr = dr / 10.0;
    // dr=pow(10,-1*exponent);

    // std::cout<<dr<<" ";

    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{dr, 0, 0};
    // geometry->refreshQuantities();
    // V=geometry->totalVolume();
    // D_P=-P0*(V-V_bar)/V_bar/V_bar;
    E_bead_front = Bead_1.Energy();
    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{dr, 0, 0};
    // geometry->refreshQuantities();

    E_bead_back = Bead_1.Energy();
    grad.x = (E_bead_front - E_bead_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    // V=geometry->totalVolume();
    // D_P=-P0*(V-V_bar)/V_bar/V_bar;
    E_bead_back = Bead_1.Energy();
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, dr, 0};
    // geometry->refreshQuantities();
    // V=geometry->totalVolume();
    // D_P=-P0*(V-V_bar)/V_bar/V_bar;
    E_bead_front = Bead_1.Energy();
    grad.y = (E_bead_front - E_bead_back) / (2 * dr);

    geometry->inputVertexPositions[index] = initial_pos[index] - Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    E_bead_back = Bead_1.Energy();
    geometry->inputVertexPositions[index] = initial_pos[index] + Vector3{0, 0, dr};
    // geometry->refreshQuantities();
    // V=geometry->totalVolume();
    // D_P=-P0*(V-V_bar)/V_bar/V_bar;
    E_bead_front = Bead_1.Energy();
    grad.z = (E_bead_front - E_bead_back) / (2 * dr);

    Gradients[exponent] = grad;
    Gradient_file << Gradients[exponent].x << " " << Gradients[exponent].y << " " << Gradients[exponent].z << " \n";
  }
}

EdgeData<double> Mem3DG::Edge_sizing(VertexData<double> Vert_sizings)
{

  EdgeData<double> sizing(*mesh, 0.0);

  Vertex v1;
  Vertex v2;
  for (Edge e : mesh->edges())
  {
    v1 = e.halfedge().vertex();
    v2 = e.halfedge().next().vertex();
    sizing[e] = geometry->edgeLength(e) * std::sqrt((Vert_sizings[v1] + Vert_sizings[v2]) / 2.0);
  }
  return sizing;
}

void Mem3DG::Save_mesh(size_t current_t)
{

  // Build member variables: mesh, geometry
  Vector3 Pos;

  // Here is where we need to find the latest directory
  std::string test_name = basic_name + "Ip_mem_" + std::to_string(current_t) + ".obj";
  const char *dir = test_name.c_str();
  struct stat sb;
  size_t iteration = current_t;
  int status = -1;
  while (status < 0)
  {

    if (stat(dir, &sb) != 0)
    {
      // std::cout<<"Path already exists
      status = -1;
      break;
    }
    iteration += 1;
    test_name = basic_name + "Ip_mem_" + std::to_string(iteration) + ".obj";
    dir = test_name.c_str();
  }
  std::ofstream o(test_name);
  o << "#This is a meshfile from a saved state\n";

  for (Vertex v : mesh->vertices())
  {
    Pos = geometry->inputVertexPositions[v];
    o << "v " << Pos.x << " " << Pos.y << " " << Pos.z << "\n";
  }

  // I need to save the faces now

  for (Face f : mesh->faces())
  {
    o << "f";

    for (Vertex v : f.adjacentVertices())
    {
      o << " " << v.getIndex() + 1;
    }
    o << "\n";
  }

  // I need to save the bead positions too;

  for (size_t b = 0; b < Beads.size(); b++)
  {

    test_name = basic_name + "Ip_bead_" + std::to_string(b) + "_data.txt";
    std::ofstream bead_file(test_name, std::ios_base::app);
    bead_file << Beads[b]->Pos.x << " " << Beads[b]->Pos.y << " " << Beads[b]->Pos.z << "\n";
    bead_file.close();
  }

  return;
}

int Mem3DG::remesh(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, RemeshOptions options)
{
  geometrycentral::surface::MutationManager mm(mesh, geom);
  return this->remesh(mesh, geom, mm, options);
}

int Mem3DG::remesh(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm, RemeshOptions options)
{
  // std::cout << "Here\n";
  geom.requireFaceSizing();
  for (Face f : mesh.faces())
  {
    geom.faceSizing[f] = clamp(geom.faceSizing[f] / (options.refine_angle * options.refine_angle),
                               1.0 / (options.max_absolute_length * options.max_absolute_length),
                               1.0 / (options.min_absolute_length * options.min_absolute_length));
  }

  geom.requireVertexSizing();

  bool doConnectivityChanges = true;

  options.maxIterations = 1;
  for (size_t iIt = 0; iIt < options.maxIterations; iIt++)
  {
    size_t nFlips = 10;
    if (doConnectivityChanges)
    {

      splitWorstEdges(mesh, geom, mm, options);
      nFlips = fixDelaunay(mesh, geom, mm);
      options.numberOp += nFlips;
      improveFaces(mesh, geom, mm, options);
    }

    nFlips = fixDelaunay(mesh, geom, mm);
    // std::cout<<"The number of flips is "<< nFlips<<"\n";
    options.numberOp += nFlips;

    geom.inputVertexPositions = geom.vertexPositions;
    mesh.compress();
    geom.refreshQuantities();
    double smoothing = smoothByLaplacian(mesh, geom, mm);
    while (smoothing > 1e-4)
      smoothing = smoothByLaplacian(mesh, geom, mm);
  }
  geom.unrequireFaceSizing();
  geom.unrequireVertexSizing();
  geom.inputVertexPositions = geom.vertexPositions;
  geom.purgeQuantities();

  return options.numberOp;
}

Vector3 vertexNormal(VertexPositionGeometry &geom, Vertex v, MutationManager &mm)
{
  Vector3 totalNormal = Vector3::zero();
  for (Corner c : v.adjacentCorners())
  {
    Vector3 cornerNormal = geom.cornerAngle(c) * geom.faceNormal(c.face());
    totalNormal += cornerNormal;
  }
  return normalize(totalNormal);
}

Vector3 boundaryVertexTangent(VertexPositionGeometry &geom, Vertex v, MutationManager &mm)
{
  if (v.isBoundary())
  {
    auto edgeVec = [&](Edge e) -> Vector3
    {
      return (geom.vertexPositions[e.halfedge().tipVertex()] - geom.vertexPositions[e.halfedge().tailVertex()])
          .normalize();
    };

    Vector3 totalTangent = Vector3::zero();
    for (Edge e : v.adjacentEdges())
    {
      if (e.isBoundary())
      {
        totalTangent += edgeVec(e);
      }
    }
    return totalTangent.normalize();
  }
  else
  {
    return Vector3::zero();
  }
}

inline Vector3 projectToPlane(Vector3 v, Vector3 norm) { return v - norm * dot(norm, v); }

inline Vector3 projectToLine(Vector3 v, Vector3 tangent) { return tangent * dot(tangent, v); }

bool isDelaunay_improv(VertexPositionGeometry &geom, Edge e)
{
  Halfedge he = e.halfedge();
  Vector3 p0 = geom.vertexPositions[he.vertex()];
  Vector3 p1 = geom.vertexPositions[he.next().vertex()];
  Vector3 p2 = geom.vertexPositions[he.next().next().vertex()];
  Vector3 p3 = geom.vertexPositions[he.twin().next().next().vertex()];

  double la = (p0 - p1).norm();
  double lb = (p1 - p2).norm();
  double lc = (p2 - p0).norm();
  double ld = (p3 - p1).norm();
  double le = (p3 - p0).norm();

  return (lb * lb + lc * lc - la * la) / (lb * lc) + (ld * ld + le * le - la * la) / (ld * le) >= -0.1;
}

inline double diamondAngle(Vector3 a, Vector3 b, Vector3 c, Vector3 d) // dihedral angle at edge a-b
{
  Vector3 n1 = cross(b - a, c - a);
  Vector3 n2 = cross(b - d, a - d);
  return PI - angle(n1, n2);
}

inline bool checkFoldover(Vector3 a, Vector3 b, Vector3 c, Vector3 x, double angle)
{
  return diamondAngle(a, b, c, x) < angle;
}

inline Vector3 edgeMidpoint(SurfaceMesh &mesh, VertexPositionGeometry &geom, Edge e)
{
  Vector3 endPos1 = geom.vertexPositions[e.halfedge().tailVertex()];
  Vector3 endPos2 = geom.vertexPositions[e.halfedge().tipVertex()];
  return (endPos1 + endPos2) / 2;
}

Vector3 findCircumcenter(Vector3 p1, Vector3 p2, Vector3 p3)
{
  // barycentric coordinates of circumcenter
  double a = (p3 - p2).norm();
  double b = (p3 - p1).norm();
  double c = (p2 - p1).norm();
  double a2 = a * a;
  double b2 = b * b;
  double c2 = c * c;
  Vector3 circumcenterLoc{a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2)};
  // normalize to sum of 1
  circumcenterLoc = normalizeBarycentric(circumcenterLoc);

  // change back to space
  return circumcenterLoc[0] * p1 + circumcenterLoc[1] * p2 + circumcenterLoc[2] * p3;
}

Vector3 findCircumcenter(VertexPositionGeometry &geom, Face f)
{
  // retrieve the face's vertices
  int index = 0;
  Vector3 p[3];
  for (Vertex v0 : f.adjacentVertices())
  {
    p[index] = geom.vertexPositions[v0];
    index++;
  }
  return findCircumcenter(p[0], p[1], p[2]);
}

Vector3 findODTCenter(VertexPositionGeometry &geom, Face f, MutationManager &mm)
{
  Vector3 p0 = geom.vertexPositions[f.halfedge().tailVertex()];
  Vector3 p1 = geom.vertexPositions[f.halfedge().tipVertex()];
  Vector3 p2 = geom.vertexPositions[f.halfedge().next().tipVertex()];

  for (Edge e : f.adjacentEdges())
  {
    if (e.isBoundary() || !mm.mayFlipEdge(e))
    {
      // e is not flippable. return barycenter
      return (p0 + p1 + p2) / 3.;
    }
  }
  return findCircumcenter(p0, p1, p2);
}

bool shouldCollapse(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, Edge e)
{
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges())
  {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2)
    {
      edgesToCheck.push_back(he.next());
    }
  }

  // link around v2
  for (Halfedge he : v2.outgoingHalfedges())
  {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1)
    {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);
  // Vector3 butterfly = edgeButterfly(mesh,geom,e);
  for (Halfedge he0 : edgesToCheck)
  {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2))
    {
      return false;
    }
  }

  return true;
}

bool shouldCollapse(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, Edge e, RemeshOptions options)
{
  std::vector<Halfedge> edgesToCheck;
  Vertex v1 = e.halfedge().vertex();
  Vertex v2 = e.halfedge().twin().vertex();

  // find (halfedge) link around the edge, starting with those surrounding v1
  for (Halfedge he : v1.outgoingHalfedges())
  {
    if (he.next().tailVertex() != v2 && he.next().tipVertex() != v2)
    {
      edgesToCheck.push_back(he.next());
    }
  }

  // link around v2
  for (Halfedge he : v2.outgoingHalfedges())
  {
    if (he.next().tailVertex() != v1 && he.next().tipVertex() != v1)
    {
      edgesToCheck.push_back(he.next());
    }
  }

  // see if the point that would form after a collapse would cause a major foldover with surrounding edges
  double area;
  double perimeter;
  double aspect;
  double a0;
  Vector3 midpoint = edgeMidpoint(mesh, geom, e);
  double E_metric;
  double NewSizing = std::max(geom.vertexSizing[v1], geom.vertexSizing[v2]);
  // Vector3 butterfly = edgeButterfly(mesh,geom,e);
  for (Halfedge he0 : edgesToCheck)
  {
    Halfedge heT = he0.twin();
    Vertex v1 = heT.tailVertex();
    Vertex v2 = heT.tipVertex();
    Vertex v3 = heT.next().tipVertex();
    Vector3 a = geom.vertexPositions[v1];
    Vector3 b = geom.vertexPositions[v2];
    Vector3 c = geom.vertexPositions[v3];
    if (checkFoldover(a, b, c, midpoint, 2))
    {
      return false;
    }

    v1 = he0.next().tipVertex();
    v2 = he0.tailVertex();
    v3 = he0.tipVertex();

    a = midpoint;
    b = geom.vertexPositions[v2];
    c = geom.vertexPositions[v3];

    a0 = geom.faceAreas[he0.face()];
    area = 0.5 * norm(cross(b - a, c - a));
    perimeter = norm(b - a) + norm(c - a) + norm(c - b);
    aspect = 12 * sqrt(3) * area / (perimeter * perimeter);
    if ((area < a0 && area < 0.1 * options.min_absolute_length * options.min_absolute_length) ||
        aspect < options.aspect_min)
      return false;

    // We check the metric is not too big
    heT = he0;
    if (geom.edgeLengths[heT.edge()] < 1e-10 || geom.edgeLengths[heT.next().edge()] < 1e-10 ||
        geom.edgeLengths[heT.next().next().edge()] < 1e-10)
      return false;

    // Ok we need to change this

    heT = heT.next();

    E_metric = norm(c - a) * sqrt((NewSizing + geom.vertexSizing[heT.tailVertex()]) / 2.0);
    if (E_metric > 0.9)
      return false;

    heT = heT.next();
    E_metric = norm(b - a) * sqrt((NewSizing + geom.vertexSizing[heT.tipVertex()]) / 2.0);
    if (E_metric > 0.9)
      return false;

    heT = heT.next();
    E_metric = norm(b - c) * sqrt((geom.vertexSizing[heT.tailVertex()] + geom.vertexSizing[heT.tipVertex()]) / 2.0);
    if (E_metric > 0.9)
      return false;
  }

  return true;
}

struct Deterministic_sort2
{
  inline bool operator()(const std::pair<double, Edge> &left, const std::pair<double, Edge> &right)
  {
    return left.first < right.first;
  }
} deterministic_sort2;

std::vector<Edge> findBadEdges(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                               RemeshOptions options)
{

  std::vector<std::pair<double, Edge>> edgems;
  for (Edge e : mesh.edges())
  {
    // TODO add the no remesh option here (not necessary yet) This is the spirit
    // if (options.no_remesh_list) {
    //   for (Edge e : mesh.edges()) {

    //     if (options.No_remesh_list[e] == 0) {

    //       toSplit.push_back(e);
    //     }
    //   }
    // } else {
    //   for (Edge e : mesh.edges()) {
    //     toSplit.push_back(e);
    //   }
    // }
    double E_sizing = geom.edgeLengths[e] *
                      sqrt(geom.vertexSizing[e.halfedge().vertex()] + geom.vertexSizing[e.halfedge().twin().vertex()]) /
                      (2.0);
    if (E_sizing > 1 && geom.edgeLengths[e] > 0.06)
      edgems.push_back(std::make_pair(E_sizing, e));
  }

  std::sort(edgems.begin(), edgems.end(), deterministic_sort2);
  std::vector<Edge> edges(edgems.size());
  for (size_t e = 0; e < edgems.size(); e++)
  {
    edges[e] = edgems[edgems.size() - e - 1].second;
  }
  return edges;
}

bool Mem3DG::splitWorstEdges(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                             RemeshOptions options)
{
  geom.requireVertexDualAreas();
  geom.requireVertexSizing();
  geom.requireEdgeLengths();

  bool didSplit = false;
  std::vector<Edge> toSplit = findBadEdges(mesh, geom, mm, options);
  std::vector<Face> activeFaces;

  // Here we will check
  for (Edge e : toSplit)
  {
    if (geom.edgeLengths[e] < 0.06)
    {
      std::cout << "THere is an edge with unacceptable size in the queue\n";
    }
  }

  double newSizing;
  while (!toSplit.empty())
  {
    Edge e = toSplit.back();
    toSplit.pop_back();
    double length_e = geom.edgeLength(e);

    newSizing = (geom.vertexSizing[e.halfedge().vertex()] + geom.vertexSizing[e.halfedge().twin().vertex()]) / (2.0);

    Vector3 newPos = edgeMidpoint(mesh, geom, e);
    Halfedge he = mm.splitEdge(e, newPos);
    if (he != Halfedge())
    {
      geom.vertexSizing[he.vertex()] = newSizing;
      Halfedge heround = he;
      int counter = 0;
      do
      {
        activeFaces.push_back(heround.face());
        counter += 1;
        if (geom.edgeLengths[heround.edge()] < 1e-10)
        {
          geom.edgeLengths[heround.edge()] =
              norm(geom.vertexPositions[heround.vertex()] - geom.vertexPositions[heround.tipVertex()]);
        }
        heround = heround.twin().next();
      } while (heround != he);

      didSplit = true;
      options.numberOp += 1;
      // flipSubset(activeFaces, mesh, geom, mm, options);
    }
  }

  // We splitted everything
  // return;
  // std::cout<<"THe number of edges that i shall not collapse are"<<counter_vertex<<" this number shouuld be
  // constant?\n"; actually collapsing
  geom.unrequireEdgeLengths();
  geom.unrequireVertexDualAreas();
  geom.unrequireVertexSizing();
  mesh.compress();
  return didSplit;
}

bool Mem3DG::improveFaces(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm, RemeshOptions options)
{
  geom.requireVertexDualAreas();
  geom.requireVertexSizing();
  geom.requireEdgeLengths();

  bool didCollapse = false;
  // queues of edges to CHECK to change
  std::vector<Edge> toCollapse;

  if (options.no_remesh_list)
  {
    for (Edge e : mesh.edges())
    {

      if (options.No_remesh_list[e] == 0)
      {
        // std::cout<<"Not remeshing this edge "<<e.getIndex() <<"\n";
        toCollapse.push_back(e);
      }
    }
  }
  else
  {
    for (Edge e : mesh.edges())
    {
      toCollapse.push_back(e);
    }
  }

  size_t counter_vertex = 0;
  // actually splitting
  double E_sizing;
  double newSizing;

  while (!toCollapse.empty())
  {

    Edge e = toCollapse.back();
    toCollapse.pop_back();
    if (e == Edge() || e.isDead())
      continue; // make sure it exists

    // Now we do
    if (geom.edgeLengths[e] < options.max_absolute_length)
    {
      if (geom.edgeLengths[e] < 1e-4)
        std::cout << "THe edgelength check says " << geom.edgeLengths[e] << " \n";
      Vector3 newPos = edgeMidpoint(mesh, geom, e);
      newSizing = std::max(geom.vertexSizing[e.halfedge().tipVertex()], geom.vertexSizing[e.halfedge().tailVertex()]);
      if (shouldCollapse(mesh, geom, e, options))
      {
        Vertex v = mm.collapseEdge(e, newPos);
        if (v != Vertex())
        {
          options.numberOp += 1;
          geom.vertexSizing[v] = newSizing;
          // std::vector<Face> active_faces;
          // for (Face f : v.adjacentFaces()) active_faces.push_back(f);
          // flipSubset(active_faces, mesh, geom, mm, options);
          didCollapse = true;
        }
      }
    }
  }
  geom.unrequireEdgeLengths();
  geom.unrequireVertexDualAreas();
  geom.unrequireVertexSizing();
  mesh.compress();
  return didCollapse;
}

/*
  This function is not actually working x.x the Energy criteria is coded but doesnt work :p
*/
size_t Mem3DG::fixDelaunay(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm)
{
  // Logic duplicated from surface/intrinsic_triangulation.cpp

  // return 0;
  std::deque<Edge> edgesToCheck;      // queue of edges to check if Delaunay
  EdgeData<bool> inQueue(mesh, true); // true if edge is currently in edgesToCheck
  Halfedge he;
  Halfedge heS;
  bool VConstraint = false;
  int VConstraintIdx = 0;
  double VolContrib = 0.0;
  double NewVolContrib = 0.0;
  // start with all edges
  for (Edge e : mesh.edges())
  {
    edgesToCheck.push_back(e);
  }

  // First thing i need to do is get the total energy right. Or better yet.
  std::vector<double> vertexEnergies(mesh.nVertices(), 0.0);
  std::vector<double> Ones(mesh.nVertices(), 1.0);

  // Function that calculates the energy per vertex
  int BCounter = 0;
  VertexData<double> vertexE(mesh, 0.0);
  for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
  {
    if (Sim_handler->Energies[i] == "Bending")
      vertexE += Sim_handler->Ev_Bending(Sim_handler->Energy_constants[i]);
    if (Sim_handler->Energies[i] == "Bending_tan")
      vertexE += Sim_handler->Ev_Bending_tan(Sim_handler->Energy_constants[i]);
    if (Sim_handler->Energies[i] == "Surface_Tension")
      vertexE += Sim_handler->Ev_SurfaceTension(Sim_handler->Energy_constants[i]);
    if (Sim_handler->Energies[i] == "Volume_Constraint")
    {
      VConstraint = true;
      VConstraintIdx = i;
    }
    if (Sim_handler->Energies[i] == "Bead")
    {
      // This is gonna be a hustle
      vertexE += Sim_handler->Beads[BCounter]->Bead_I->V_Tot_Energy();
      BCounter += 1;
    }
  }
  BCounter = 0;
  // Ok so these are the vertex energies
  double Current_vol = geometry->totalVolume();
  double Current_area = geometry->totalArea();

  bool delaunay = true;
  // OK great next thing is to calculate the energy. So Lets see
  // Bending energy can be calculated per vertex easy.
  // Surface tension can be calculated per vertex using dual areas
  // Volume and area constraints can only be calculated as totals.

  // counter and limit for number of flips
  size_t flipMax = 100 * mesh.nVertices();
  size_t nFlips = 0;

  VertexData<int> CheckedV(mesh, 0);
  std::vector<Vertex> VertToCheck(0);
  std::vector<double> VertNewE(0);
  // bool checkNeigh = false;
  // return 0;
  while (!edgesToCheck.empty() && nFlips < flipMax)
  {
    VertToCheck.resize(0);
    VertNewE.resize(0);
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    if (e.isBoundary())
      continue;

    // Sooo i need to decide which vertices will get explored
    he = e.halfedge();

    // VertToCheck.push_back(he.vertex());

    for (int s = 0; s < 2; s++)
    {

      heS = he.next();
      for (int f = 0; f < 2; f++)
      {
        for (Vertex v : heS.twin().face().adjacentVertices())
        {
          if (CheckedV[v] == 0)
          {
            VertToCheck.push_back(v);
            VertNewE.push_back(0.0);
            CheckedV[v] = 1;
          }
        }
        heS = heS.next();
      }
      he = he.twin();
    }

    // So now i have all the vertices i want to check

    // Next thing is to get the current energy of this configuration, but actually we already have that
    double Curr_E = 0.0;
    double New_E = 0.0;
    for (Vertex v : VertToCheck)
    {
      Curr_E += vertexE[v];
    }
    // Ok and the volume if the constraint allows it
    if (VConstraint)
    {
      VolContrib = 0.0;
      // I need to put the volume contribution
      he = e.halfedge();
      VolContrib = geometry->faceVolume(he.face()) + geometry->faceVolume(he.twin().face());
      //  In this case both faces contribute to the volume.
    }

    // Then we flip

    bool ECheckFLip = mm.flipEdge(e);
    // Ok its time to check the energy of the new configuration.
    if (ECheckFLip)
    {
      BCounter = 0;
      for (size_t i = 0; i < Sim_handler->Energies.size(); i++)
      {
        if (Sim_handler->Energies[i] == "Bending")
        {
          for (size_t j = 0; j < VertToCheck.size(); j++)
          {
            Vertex v = VertToCheck[j];
            double val = Sim_handler->V_Bending(Sim_handler->Energy_constants[i], v);
            VertNewE[j] += val;
            New_E += val;
          }
        }
        if (Sim_handler->Energies[i] == "Bending_tan")
        {
          for (size_t j = 0; j < VertToCheck.size(); j++)
          {
            Vertex v = VertToCheck[j];
            double val = Sim_handler->V_Bending_tan(Sim_handler->Energy_constants[i], v);
            VertNewE[j] += val;
            New_E += val;
          }
        }
        if (Sim_handler->Energies[i] == "Surface_Tension")
          for (size_t j = 0; j < VertToCheck.size(); j++)
          {
            Vertex v = VertToCheck[j];
            double val = geometry->vertexDualArea(v) * Sim_handler->Energy_constants[i][0];
            VertNewE[j] += val;
            New_E += val;
          }
        if (Sim_handler->Energies[i] == "Bead")
        {
          for (size_t j = 0; j < VertToCheck.size(); j++)
          {
            Vertex v = VertToCheck[j];
            double val = Sim_handler->Beads[BCounter]->Bead_I->V_Energy(v);
            VertNewE[j] += val;
            New_E += val;
          }
          BCounter += 1;
        }
      }
      if (VConstraint)
      {
        NewVolContrib = 0.0;
        he = e.halfedge();
        NewVolContrib = geometry->faceVolume(he.face()) + geometry->faceVolume(he.twin().face());
        // We also add the volume change
        double KV = Sim_handler->Energy_constants[VConstraintIdx][0];
        double V_bar = Sim_handler->Energy_constants[VConstraintIdx][1];
      }

      // Ok
      // So i have CurrE and NewE
      if (Curr_E - New_E > 1e-3)
      {
        // In this case the flip minimizes energy
        nFlips++;

        // Now i need to rewrite the edge energies
        for (size_t k = 0; k < VertToCheck.size(); k++)
        {
          vertexE[VertToCheck[k]] = VertNewE[k];
          CheckedV[VertToCheck[k]] = 0;
          Current_vol = Current_vol - VolContrib + NewVolContrib;

          // Here we will reset everything
        }
      }
      else
      {

        // If the change in energy is not small and the
        // We also want to keep this one
        mm.flipEdge(e);
        continue;
      }
    }
    else
    {
      continue;
    }

    // bool delaunay = isDelaunay_improv(geom, e);

    // if not Delaunay, try to flip edge
    // bool wasFlipped = mm.flipEdge(e);

    // if (!wasFlipped)
    //   continue;

    // nFlips++;

    // Add neighbors to queue, as they may need flipping now
    // he = e.halfedge();
    // std::array<Edge, 4> neighboringEdges{he.next().edge(), he.next().next().edge(), he.twin().next().edge(),
    //                                      he.twin().next().next().edge()};
    // for (Edge nE : neighboringEdges)
    // {
    //   if (!inQueue[nE])
    //   {
    //     edgesToCheck.push_back(nE);
    //     inQueue[nE] = true;
    //   }
    // }
  }
  return nFlips;
}

size_t fixDelaunay(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                   RemeshOptions options)
{
  // Logic duplicated from surface/intrinsic_triangulation.cpp

  std::deque<Edge> edgesToCheck;      // queue of edges to check if Delaunay
  EdgeData<bool> inQueue(mesh, true); // true if edge is currently in edgesToCheck

  // start with all edges
  for (Edge e : mesh.edges())
  {
    edgesToCheck.push_back(e);
  }
  // We check all the vertices.
  VertexData<double> vertexE(mesh);

  double E = 0;
  // Calculate the energy per vertex

  // counter and limit for number of flips
  size_t flipMax = 100 * mesh.nVertices();
  size_t nFlips = 0;
  while (!edgesToCheck.empty() && nFlips < flipMax)
  {
    Edge e = edgesToCheck.front();
    edgesToCheck.pop_front();
    inQueue[e] = false;

    // Now we do a function that checks if the ENergy decreases
    if (e.isBoundary() || isDelaunay_improv(geom, e))
      continue;

    // if not Delaunay, try to flip edge
    bool wasFlipped = mm.flipEdge(e);

    if (!wasFlipped)
      continue;

    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge he = e.halfedge();
    std::array<Edge, 4> neighboringEdges{he.next().edge(), he.next().next().edge(), he.twin().next().edge(),
                                         he.twin().next().next().edge()};
    for (Edge nE : neighboringEdges)
    {
      if (!inQueue[nE])
      {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }
  return nFlips;
}

double Mem3DG::smoothByLaplacian(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm, double stepSize,
                                 RemeshBoundaryCondition2 bc)
{
  VertexData<Vector3> vertexOffsets(mesh);
  geom.requireVertexNormals();
  geom.requireVertexPositions();

  for (Vertex v : mesh.vertices())
  {
    // calculate average of surrounding vertices
    Vector3 avgNeighbor = Vector3::zero();
    for (Vertex j : v.adjacentVertices())
    {
      avgNeighbor += geom.vertexPositions[j];
    }
    avgNeighbor /= v.degree();

    Vector3 updateDirection = avgNeighbor - geom.vertexPositions[v];

    // project updateDirection onto space of allowed movements
    Vector3 stepDir;
    if (v.isBoundary())
    {
      switch (bc)
      {
      case RemeshBoundaryCondition2::Fixed:
        stepDir = Vector3::zero();
        break;
      case RemeshBoundaryCondition2::Tangential:
        // for free boundary vertices, project the average to the boundary tangent line
        stepDir = projectToLine(updateDirection, boundaryVertexTangent(geom, v, mm));
        break;
      case RemeshBoundaryCondition2::Free:
        // for free boundary vertices, project the average to the surface tangent plane
        stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
        break;
      }
    }
    else
    {
      // for interior vertices, project the average to the tangent plane
      stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
    }
    vertexOffsets[v] = stepSize * stepDir;
  }

  // update final vertices
  double totalMovement = 0;
  for (Vertex v : mesh.vertices())
  {
    bool didMove = mm.repositionVertex(v, vertexOffsets[v]);
    if (didMove)
    {
      totalMovement += vertexOffsets[v].norm();
    }
  }
  // geom.inputVertexPositions = geom.vertexPositions;
  geom.unrequireVertexNormals();
  geom.unrequireVertexPositions();
  // What i would like to check is if there is any nAN VERTEX AFTERTHIS
  return totalMovement / mesh.nVertices();
}

double Mem3DG::smoothByCircumcenter(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geom, MutationManager &mm,
                                    double stepSize, RemeshBoundaryCondition2 bc)
{
  geom.requireFaceAreas();
  VertexData<Vector3> vertexOffsets(mesh);
  for (Vertex v : mesh.vertices())
  {
    Vector3 updateDirection = Vector3::zero();
    for (Face f : v.adjacentFaces())
    {
      // add the circumcenter weighted by face area to the update direction
      Vector3 circum = findODTCenter(geom, f, mm);
      updateDirection += geom.faceArea(f) * (circum - geom.vertexPositions[v]);
    }
    updateDirection /= (6 * geom.vertexDualArea(v));

    // project updateDirection onto space of allowed movements
    Vector3 stepDir;
    if (v.isBoundary())
    {
      switch (bc)
      {
      case RemeshBoundaryCondition2::Fixed:
        stepDir = Vector3::zero();
        break;
      case RemeshBoundaryCondition2::Tangential:
        // for free boundary vertices, project the average to the boundary tangent line
        stepDir = projectToLine(updateDirection, boundaryVertexTangent(geom, v, mm));
        break;
      case RemeshBoundaryCondition2::Free:
        // for free boundary vertices, project the average to the surface tangent plane
        stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
        break;
      }
    }
    else
    {
      // for interior vertices, project the average to the tangent plane
      stepDir = projectToPlane(updateDirection, vertexNormal(geom, v, mm));
    }
    vertexOffsets[v] = stepSize * stepDir;
  }

  // update final vertices
  double totalMovement = 0;
  for (Vertex v : mesh.vertices())
  {
    bool didMove = mm.repositionVertex(v, vertexOffsets[v]);
    if (didMove)
    {
      totalMovement += vertexOffsets[v].norm();
    }
  }
  return totalMovement / mesh.nVertices();
}

// ---- from Mem-3dg.cpp, around line 1050
  // for (size_t bi = 0; bi < Beads.size(); bi++)
  // {
  //   if (Beads[bi]->state == "manual")
  //   {
  //     Vector3 Bpos = Beads[bi]->Pos;
  //     if ((Bpos.norm2() < 4.0 && dot(Bpos, Beads[bi]->Velocity) < 0) || (Bpos.norm2() > 4.0 && dot(Bpos, Beads[bi]->Velocity) > 0)) // The 2.0 here is hardcoded and it means the radius of the vesicle
  //     {
  //       std::cout << "\t\t Manual bead because it moved too much\n";
  //       std::cout << "The bead positions 2 is" << sqrt(Bpos.norm2()) << " \n";
  //       Beads[bi]->state = "default";
  //     }
  //   }
  // }

// ---- from Mem-3dg.cpp, around line 2198
      // I can do this i have grad and grad_theory so i can actually compare them

      // I want to know a little more abt this direction.

      // r= Bead_1.Pos- geometry->inputVertexPositions[v];
      // r_dist=r.norm();
      // r= r.unit();
      // Area_grad=Grad_area[v].unit();
      // // Area_grad= geometry->vertexNormalMeanCurvature(v).unit();

      // double E_v=4*1.0*(pow(1.0/r_dist,12)-pow(1.0/r_dist,6));
      // Vector3 F2=E_v *-1*geometry->vertexNormalMeanCurvature(v);

