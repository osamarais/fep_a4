// This file contains all the functions for the element routines

// Header
std::pair<bool,double> essential_BC(int boundary_number);
std::pair<double,double> natural_BC(int boundary_number);


struct contribution{
  int row;
  int column;
  double coefficient;
  double known;
};

struct boundary_struct{
  int boundary;
  pMeshEnt e;
};



void const_str_tri(pMeshEnt e);

// Minimum requirements are 3 sided and 4 sided Linear and Quaratic Elements

// The type of elements in the mesh must be identified. Functions for that are placed here
/*

Element enumeration is as follows:
1: Linear Tet
2: Quadratic (Serendipity) Tet
3: Linear Quad
4: Quadratic (Serendipity) Quad

*/

int get_element_type(pMesh mesh, pMeshEnt e){
  int element_type = 0;
  int edgenodes = 0;
  int numberofedges = 0;
  // get nummber of nodes on the edges
  pShape myshape = pumi_mesh_getShape(mesh);
  edgenodes = pumi_shape_getNumNode(myshape, 1);
  // get number of nodes on the faces
  // implement this later if needed
  // check the nmber of adjacent edges to get the shape of the element
  numberofedges = pumi_ment_getNumAdj(e,1);
  // Now use the data about the number of edge nodes etc. to get the element type
  if ((0 == edgenodes) && (3 == numberofedges)){
    element_type = 1;
  }
  else if ((1 == edgenodes) && (3 == numberofedges)){
    element_type = 2;
  }
  else if ((0 == edgenodes) && (4 == numberofedges)){
    element_type = 3;
  }
  else if ((1 == edgenodes) && (4 == numberofedges)){
    element_type = 4;
  }
  else {
    element_type = 0;
  }

  if (0 == element_type){
    printf("Unidentified Element!!!\n");
  }
  //printf("The element_type type identified is %d \n", element_type);
  return element_type;
}








// Get all the mesh edges classified on ALL geometric boundaries


void get_all_boundary_edges(pGeom &g,pMesh &mesh, std::vector<boundary_struct> &mesh_ents){
  // Iterate over all the geometric edges
  // put mesh entities classified on those edges into the container
  for (pGeomIter it1 = g->begin(1); it1!=g->end(1);++it1){
    pGeomEnt ge = *it1;
    std::vector<pMeshEnt> entities_on_this;
    // get all the reverse classified entites on this
    pumi_gent_getRevClas(ge, entities_on_this);
    // push the subcontainer into the main container
    for (std::vector<pMeshEnt>::iterator it2 = entities_on_this.begin(); it2!= entities_on_this.end(); ++it2){
      int faceid = pumi_gent_getID(*it1);
      boundary_struct this_edge;
      this_edge.boundary = pumi_gent_getID(*it1);
      this_edge.e = *it2;
      mesh_ents.push_back(this_edge);
    }
  }
}


// Get a list of all the vertices on boundaries

void get_all_boundary_nodes(pMesh &mesh, std::vector<boundary_struct> &boundary_edges, std::vector<boundary_struct> &mesh_ents){
  // Get the adjacent vertices of all the edges.
  // avoid duplication by checking the list and tallying against e and boundary both
  // allow vertex to be classifoed on more than one boundary
  // Iterate over the edges
  for (std::vector<boundary_struct>::iterator it1 = boundary_edges.begin(); it1!= boundary_edges.end(); ++it1){
    boundary_struct this_edge = *it1;
    // Get the adjacent edges
    // go through the container, check if both boundary and e have not been satisfied. if not push back
    Adjacent adjacent;
    pumi_ment_getAdjacent(this_edge.e,0,adjacent);
    if (mesh_ents.size() == 0){
      // push back the first vertices to intitialize the iteration
      boundary_struct this_vert;
      this_vert.boundary = this_edge.boundary;
      this_vert.e = adjacent[0];
      mesh_ents.push_back(this_vert);
      //printf("Added first vertex\n");
    }
    // Iterate over adjacents to check of they are in the mesh_ents or not
    for(int i = 0; i< adjacent.size(); i ++){
      bool found = false;
      boundary_struct this_vert;
      for (std::vector<boundary_struct>::iterator it2 = mesh_ents.begin(); it2!= mesh_ents.end(); ++it2){
        //printf("Finding something here\n");
        // if found flag and break
        this_vert.e = adjacent[i];
        this_vert.boundary = this_edge.boundary;
        boundary_struct other_vert;
        other_vert = *it2;
        if ((this_vert.boundary == other_vert.boundary)&&(this_vert.e == other_vert.e)){
          found = true;
          break;
        }
      }
      // if flag is still false, add the vertex
      if (!found){
        mesh_ents.push_back(this_vert);
        //printf("Vertex %d has been added \n", pumi_ment_getID(this_vert.e));
      }
    }
  }
}


// Function returns which boundaries a certain mesh entity is on

void get_bound_num(pMeshEnt ment, std::vector<boundary_struct> &mesh_ents, std::vector<int> &list){
  for(std::vector<boundary_struct>::iterator it1 = mesh_ents.begin() ; it1!=mesh_ents.end(); ++it1){
    boundary_struct this_element = *it1;
    if (this_element.e == ment){
      list.push_back(this_element.boundary);
    }
  }
}


// In this function, pass a mesh entity (face) and get the area of that face

double get_face_area(pMeshEnt face){
  // Get all the adjacent vertices
  Adjacent vertices;
  int num_verts = pumi_ment_getAdjacent(face,0,vertices);
  // Get the coordinates of the vertices
  //std::vector<double*> coordinates = std::vector<double*>();
  double area = 0;
  // create determinant structure
  for (int i=0; i<num_verts; i++){
    if(i<num_verts-1){
      double xyz1[3];
      double xyz2[3];
      pumi_node_getCoord(vertices[i],0,xyz1);
      pumi_node_getCoord(vertices[i+1],0,xyz2);
      area = area + 0.5*( (xyz1[0]*xyz2[1]) - (xyz1[1]*xyz2[0]) );
    }
    else{
      double xyz1[3];
      double xyz2[3];
      pumi_node_getCoord(vertices[i],0,xyz1);
      pumi_node_getCoord(vertices[0],0,xyz2);
      area = area + 0.5*( (xyz1[0]*xyz2[1]) - (xyz1[1]*xyz2[0]) );
    }
    //double xy[2];
    //xy[0] = xyz[0];
    //xy[1] = xyz[1];

    //coordinates.push_back(xyz);
  }
  return abs(area);
}






// Reorder Element to make it CCW (swapping Algorithm)
// Input: a vector of vertices
// Output: a vector of vertices
bool reorder_verts(Adjacent &adjacent){
  // keep looping till this counter reaches the size of the vector container
  std::vector<pMeshEnt> vertices;
  for (int i = 0; i < adjacent.size(); i++){
    vertices.push_back(adjacent[i]);
  }
  //printf("we are here \n");
  for (int i = 0; i < (vertices.size()-2) ; i++){
    // Get the coordinates
    double coord1[3];
    double coord2[3];
    double coord3[3];
    pumi_node_getCoord(vertices[0], 0, coord1);
    pumi_node_getCoord(vertices[i+1], 0, coord2);
    pumi_node_getCoord(vertices[i+2], 0, coord3);
    //printf("%d \n", i);
    /*
    printf("x vertex 1 %f\n", coord1[0]);
    printf("y vertex 1 %f\n", coord1[1]);
    printf("x vertex 2 %f\n", coord2[0]);
    printf("y vertex 2 %f\n", coord2[1]);
    printf("x vertex 3 %f\n", coord3[0]);
    printf("y vertex 3 %f\n", coord3[1]);
    */
    double shoelace = coord1[0]*coord2[1] + coord2[0]*coord3[1] + coord3[0]*coord1[1] - coord1[1]*coord2[0] - coord2[1]*coord3[0] - coord3[1]*coord1[0];
    //printf("shoelace is %f\n", shoelace);
    if (shoelace<0){
      // swap the elements and restart the counter
      printf("SWAPPED!!! \n");
      std::swap(vertices[i+1],vertices[i+2]);
      i = 0;
      return true;
    }
    // Keep checking the shoelace against three every times.
    // Upon failure, swap the last two and then start over i.e. i = 0

  }
  return false;
}






contribution region_routine(pMesh mesh, pMeshEnt e, std::vector<contribution> &region_contributions){
  //contribution region_contribution;
  switch (get_element_type(mesh,e)) {
    case 0:
    printf("      Unrecognised element!!!!!!!!!!!\n");
    break;
    case 1:
    printf("Constant Strain Triangle Element\n");
    //const_str_tri(e, region_contributions);
    break;
    default:
    printf("      Unknown Error in Element routine!!!!!!!!!!\n");
    break;
  }
}



void const_str_tri(pMeshEnt e, std::vector<contribution> &region_contributions){
  int id = pumi_ment_getID(e);
  printf("\n\n\n\nRegion %d:\n\n",id);
  // Get global coordinates
  Adjacent adjacent;
  pumi_ment_getAdjacent(e,0,adjacent);
  double coord1[3];
  double coord2[3];
  double coord3[3];
  pumi_node_getCoord(adjacent[0], 0, coord1);
  pumi_node_getCoord(adjacent[1], 0, coord2);
  pumi_node_getCoord(adjacent[2], 0, coord3);
  printf("x vertex 1 %f\n", coord1[0]);
  printf("y vertex 1 %f\n", coord1[1]);
  printf("x vertex 2 %f\n", coord2[0]);
  printf("y vertex 2 %f\n", coord2[1]);
  printf("x vertex 3 %f\n", coord3[0]);
  printf("y vertex 3 %f\n\n", coord3[1]);
  // Generate Jacobian
  double J[2][2];
  J[0][0] = coord1[0] - coord3[0];
  J[0][1] = coord1[1] - coord3[1];
  J[1][0] = coord2[0] - coord3[0];
  J[1][1] = coord2[1] - coord3[1];
  // Multiply and get jj
  double JJ[2][2] = {0};
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 2; k++){
        JJ[i][j] += J[i][k]*J[k][j];
      }
    }
  }
  // Create Del Matrices and Multiply
  double del[3][2] = {{1,0},{0,1},{-1,-1}};
  double delT[2][3] = {{1,0,-1},{0,1,-1}};
  double delJJ[3][2] = {0};
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 2; k++){
        delJJ[i][j] += del[i][k]*JJ[k][j];
        printf("JJ i %d j %d    %f \n",i,j,delJJ[i][j]);
      }
    }
  }
  // Create complete matrix, and also do integral (divide by 2)
  double delJJdelT[3][3] = {0};
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      for (int k = 0; k < 2; k++){
        delJJdelT[i][j] += delJJ[i][k]*delT[k][j]/2;
      }
    }
  }
  // Now simply assemble the contribution and push it back
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      contribution c;
      c.coefficient = delJJdelT[i][j];
      c.known = 0;
      c.row = pumi_ment_getID(adjacent[i]);
      c.column = pumi_ment_getID(adjacent[j]);
      region_contributions.push_back(c);
      printf("Row %d \n", c.row);
      printf("Column%d \n", c.column);
      printf("contribution coefficient %f \n", c.coefficient);
    }
  }
}




void lin_str_tri(pMeshEnt e, std::vector<contribution> &region_contributions){
  // Get coordiantes of vertices
  int id = pumi_ment_getID(e);
  printf("\n\n\n\nRegion %d:\n\n",id);
  Adjacent adjacentv;
  pumi_ment_getAdjacent(e,0,adjacentv);
  double coord1[3] = {0};
  double coord2[3] = {0};
  double coord3[3] = {0};
  pumi_node_getCoord(adjacentv[0], 0, coord1);
  pumi_node_getCoord(adjacentv[1], 0, coord2);
  pumi_node_getCoord(adjacentv[2], 0, coord3);
  printf("x vertex 1 %f\n", coord1[0]);
  printf("y vertex 1 %f\n", coord1[1]);
  printf("x vertex 2 %f\n", coord2[0]);
  printf("y vertex 2 %f\n", coord2[1]);
  printf("x vertex 3 %f\n", coord3[0]);
  printf("y vertex 3 %f\n", coord3[1]);
  // Get coordinates of the edge nodes in proper ordering
  double coord4[3] = {0};
  double coord5[3] = {0};
  double coord6[3] = {0};
  // loop over the vertices
  // get the shared edge
  // get the node on the shared edge
  // get its vertices and put it in the coord double
  std::vector<pMeshEnt> sharededges;
  for (int i = 0; i < adjacentv.size(); i++){
    Adjacent adjacente1;
    Adjacent adjacente2;
    pumi_ment_getAdjacent(adjacentv[i],1,adjacente1);
    if (i < adjacentv.size()-1){
      pumi_ment_getAdjacent(adjacentv[i+1],1,adjacente2);
    }
    else{
      pumi_ment_getAdjacent(adjacentv[0],1,adjacente2);
    }
    //printf("Got adjacents %d\n", i);
    // Got the adjacent edges, now get the shared one
    pMeshEnt sharededge;
    for (int j = 0; j < adjacente1.size(); j++){
      for (int k = 0; k < adjacente2.size(); k++){
        //printf("Comparing adjacents  %d  %d \n", j,k);
        if (adjacente1[j]==adjacente2[k]){
          //printf("                       Found a match! \n");
          sharededge = adjacente1[j];
          sharededges.push_back(sharededge);
        }
        else {
          //printf("Failed to match \n");
        }
      }
    }
    // Now get the node on this shared edge
    // put the coords in the appropriate double
    switch (i) {
      case 0:
      pumi_node_getCoord(sharededge, 0, coord4);
      break;
      case 1:
      pumi_node_getCoord(sharededge, 0, coord5);
      break;
      case 2:
      pumi_node_getCoord(sharededge, 0, coord6);
      break;
      default:
      printf("Unknown Error!!!!\n");
    }


  }
  printf("x vertex 4 %f\n", coord4[0]);
  printf("y vertex 4 %f\n", coord4[1]);
  printf("x vertex 5 %f\n", coord5[0]);
  printf("y vertex 5 %f\n", coord5[1]);
  printf("x vertex 6 %f\n", coord6[0]);
  printf("y vertex 6 %f\n\n", coord6[1]);

  // Create x and y vectors
  double xeT[6] = {coord1[0],coord2[0],coord3[0],coord4[0],coord5[0],coord6[0]};
  double yeT[6] = {coord1[1],coord2[1],coord3[1],coord4[1],coord5[1],coord6[1]};
  // Create s and t vectors
  double s[6] = {0,1,0,0.5,0.5,0};
  double t[6] = {0,0,1,0,0.5,0.5};
  // Matrix J depends on the quadrature rule;
  // Matrices del and delT depend on the quadrature rules
  // Use 3 point triangle quadrature
  double Sin[3] = {2/3,1/6,1/6};
  double tin[3] = {1/6,2/3,1/6};
  double Win[3] = {1/3,1/3,1/3};

  // create analytical del matrix (split into the three components constant, s coefficient, and t coefficient)
  double del_analytical[6][2][3] = { {{-3,4,4},{-3,4,4}},
                          {{-1,4,0},{0,0,0}},
                          {{0,0,0},{-1,0,4}},
                          {{4,-8,-4},{0,-4,0}},
                          {{0,0,4},{0,4,0}},
                          {{0,0,-4},{4,-4,-8}} };
  // calculate this analytical del matrix using the three point quadrature
  double del[6][2] = {0};
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 3; k++){
        // Weight will not be used since it is not distributive over the matrix operations
        del[i][j] += del_analytical[i][j][0] + del_analytical[i][j][1]*Sin[k] + del_analytical[i][j][2]*tin[k];
      }
    }
  }
  // transpose del matrix
  double delT[2][6] = {0};
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 2; j++){
      delT[j][i] = del[i][j];
    }
  }
  // Create Jacobian
  double J[2][2] = {0};
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 6; k++){
        double a = 0;
        double b = 0;
        if (i == 0){
          b = del[k][0];
        }
        else {
          b = del[k][1];
        }
        if (j == 0){
          a = xeT[k];
        }
        else {
          a = yeT[k];
        }
        J[i][j] += a*b;
      }
    }
  }
  // Get Jacobian determinant
  double detJ = J[0][0]*J[1][1]-J[0][1]*J[1][0];

  // Start multiplying the matrices
  double JJ[2][2] = {0};
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 2; k++){
        JJ[i][j] += J[i][k]*J[k][j];
      }
    }
  }
  double delJJ[6][2] = {0};
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 2; k++){
        delJJ[i][j] += del[i][k]*JJ[k][j];
        //printf("delJJ i %d j %d    %f \n",i,j,delJJ[i][j]);
      }
    }
  }
  // Create complete matrix, and also do integral (divide by 2) and use weight (taking advantage of the common weight for all of them!!!!)
  // Also multiply by the Jacobian
  double delJJdelT[6][6] = {0};
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 6; j++){
      for (int k = 0; k < 2; k++){
        delJJdelT[i][j] += delJJ[i][k]*delT[k][j]/2*detJ;
      }
    }
  }

  // Now send out the assembly
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 6; j++){
      contribution c;
      c.coefficient = delJJdelT[i][j];
      c.known = 0;
      // Need to get the proper IDs of the nodes!!!
      if (i<3){
        c.row = pumi_ment_getID(adjacentv[i]);
      }
      else{
        c.row = pumi_ment_getID(sharededges[i-3]);
      }
      if (j<3){
        c.column = pumi_ment_getID(adjacentv[j]);
      }
      else{
        c.column = pumi_ment_getID(sharededges[j-3]);
      }
      //c.row = pumi_ment_getID(adjacent[i]);
      //c.column = pumi_ment_getID(adjacent[j]);
      region_contributions.push_back(c);
      printf("Row %d \n", c.row);
      printf("Column%d \n", c.column);
      printf("contribution coefficient %f \n", c.coefficient);
    }
  }

  // Quadrature is constant for the whole thing?
}












/*
// Function to get the stiffness matrix contributions of a triangular element
// The contributions will be to the vertices
// Therefore loop over the edges only
// Information needed is
// The mesh,
// The boundaries
// The boundary conditions (this is called as a function inside the function)
// return a contribution structure
struct  contribution{
int oldid;
int newid;
double coeff;
double known;
};


contribution vertex_node_subassembler(pMesh mesh, pMeshEnt e, std::vector<boundary_struct> &boundary_verts, std::vector<boundary_struct> &boundary_edges, pNumbering node_num){
// Needs: boundary lists (vertices and edges), mesh entity (node), mesh
// Returns: a contribution structure.
contribution nodecontribution;
nodecontribution.oldid = pumi_ment_getID(e);
nodecontribution.newid = pumi_node_getNumber(node_num, e, 0, 0);
// First check if a node is on a boundary.
std::vector<int> list;
get_bound_num(e, boundary_verts, list);
// Then check if an essential boundary condition has been enforced on the node at any boundary
for (int i = 0; i < list.size(); i++){
std::pair<bool,double> BC;
// Get BC
BC = essential_BC(list[i]);
if (BC.first){
nodecontribution.coeff = 1;
nodecontribution.known = BC.second;
// If so, enforce it and exit.
return nodecontribution;
}
}

// else, check the adjacent edges, and loop over them.
Adjacent adjacent;
int num_adj = pumi_ment_getAdjacent(e, 1, adjacent);
for (int i = 0; i < adjacent.size(); i++){
get_bound_num(adjacent[i], boundary_edges, list);
std::pair<double,double> BC;
if (list.size() > 0){
BC = natural_BC(list[0]);
// Get the contribution for the edge here!!
}

}


// else, check the adjacent edges, and loop over them.
// check if any of them are on a boundary.
// if any of them are on a boundary, get the contribution of that boundary.

// Now get the adjacent regions.
// Loop over the adjacent regions.
// Check the type of the region.
// Invoke the appropriate type of area contributor to get it
// Assemble all the contributions along with the renumbered number of the node.
}




contribution linear_tet(pMesh mesh, pMeshEnt e){
// Needs: a face element, mesh
// Returns: a contribution structure (id stuff ignored)

contribution tetcontribution;
// Check if this node is on a boundary where an essential boundary condition has been specified
// If on an essential boundary condition, simple enforce it.
//
// else
// get adjacent regions, loop over them, and get their contributions
//
return tetcontribution;
}


contribution lin_edge_contribuition(pMesh mesh, pMeshEnt, std::pair<double,double> BC){
contribution edge_contribution;
return edge_contribution;
}
*/











// The ordering of the elements must be ascertained properly to assemble the matrix.
// This means that: loop over the elements in any order that is found to be okay
// Simply get the numbering from the appropriate numbering we have created for the faces and the nodes.
// Over here we must check if we have already done the reordering or not so that we may use the appropriate numbering from the mesh
// Another way is to push all new element matrices into a std::vector type container (https://www.geeksforgeeks.org/vector-insert-function-in-c-stl/)
// and then sort them using some numbering. (https://www.geeksforgeeks.org/sorting-a-vector-in-c/)

// Write element routines here to get the coefficients of the elemental matrices.

// Contributor Calculators
// Stiffness Contributor
// Force Contributor (face)
// Force Contributor (edge)

// Constraint makers
// Constrain edge
// Constrain vertex


// To assemble the matrix
// 1) Loop over all the vertices
// 2) If shape is Lagrange then loop over edges and get their node too
// 3) Use the new numbering to pass these to an assembly function























//
