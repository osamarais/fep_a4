// This file contains all the functions for the element routines

// Header


struct contribution{
  int row;
  int column;
  double coefficient;
  double known;
};

struct boundary_struct{
  int boundary;
  pMeshEnt e;
  double first;     // Placeholder for either boolean or alpha
  double second;    // Placeholder for either value or beta
};

struct BC{
  double first;     // Placeholder for either boolean or alpha
  double second;    // Placeholder for either value or beta
};

BC essential_BC(int boundary_number, pMeshEnt e, pNumbering numbering);
BC natural_BC(int boundary_number);



int get_element_type(pMesh mesh, pMeshEnt e){
  // The type of elements in the mesh must be identified. Functions for that are placed here
  /*

  Element enumeration is as follows:
  1: Linear Tet
  2: Quadratic (Serendipity) Tet
  3: Linear Quad
  4: Quadratic (Serendipity) Quad

  */
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
      // fill out the natural boundary condition
      BC natural;
      natural = natural_BC(this_edge.boundary);
      this_edge.first = natural.first;
      this_edge.second = natural.second;
      // push back
      mesh_ents.push_back(this_edge);
    }
  }
}


// Get a list of ALL the nodes on boundaries
void get_all_boundary_nodes(pMesh &mesh, std::vector<boundary_struct> &boundary_edges, std::vector<boundary_struct> &mesh_ents,pNumbering numbering){
  // Get the adjacent vertices of all the edges.
  // avoid duplication by checking the list and tallying against e and boundary both
  // allow vertex to be classified on more than one boundary
  // Iterate over the edges
  for (std::vector<boundary_struct>::iterator it1 = boundary_edges.begin(); it1!= boundary_edges.end(); ++it1){
    boundary_struct this_edge = *it1;
    // Get the adjacent vertices
    // go through the container, check if both boundary and e have not been satisfied. if not push back
    Adjacent adjacent;
    pumi_ment_getAdjacent(this_edge.e,0,adjacent);
    // If edge has a node then simply push it back without checking anything since it cannot be duplicated
    if (hasNode(mesh,this_edge.e)){
      //printf("Found an edge node\n");
      // Get the BC
      boundary_struct edge_node;
      edge_node.boundary = this_edge.boundary;
      edge_node.e = this_edge.e;
      BC essential;
      essential = essential_BC(this_edge.boundary, this_edge.e, numbering);
      edge_node.first = essential.first;
      edge_node.second = essential.second;
      // push it back

      mesh_ents.push_back(edge_node);
    }


    if (mesh_ents.size() == 0){
      // push back the first vertices to intitialize the iteration
      boundary_struct this_vert;
      this_vert.boundary = this_edge.boundary;
      this_vert.e = adjacent[0];
      // Get the boundary condition
      BC essential;
      essential = essential_BC(this_vert.boundary, this_vert.e, numbering);
      this_vert.first = essential.first;
      this_vert.second = essential.second;

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
        // Get the boundary condition
        BC essential;
        essential = essential_BC(this_vert.boundary, this_vert.e,numbering);
        this_vert.first = essential.first;
        this_vert.second = essential.second;

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
bool reorder_verts(Adjacent &adjacent){
  // Input: a vector of vertices
  // Output: a vector of vertices

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


// Element routine for T3
void T3(pMeshEnt e, std::vector<contribution> &region_contributions, pNumbering numbering){
  int id = pumi_ment_getID(e);
  //printf("\n\n\n\nRegion %d:\n\n",id);
  // Get global coordinates
  Adjacent adjacent;
  pumi_ment_getAdjacent(e,0,adjacent);
  double coord1[3];
  double coord2[3];
  double coord3[3];
  pumi_node_getCoord(adjacent[0], 0, coord1);
  pumi_node_getCoord(adjacent[1], 0, coord2);
  pumi_node_getCoord(adjacent[2], 0, coord3);
  /*
  printf("x vertex 1 %f\n", coord1[0]);
  printf("y vertex 1 %f\n", coord1[1]);
  printf("x vertex 2 %f\n", coord2[0]);
  printf("y vertex 2 %f\n", coord2[1]);
  printf("x vertex 3 %f\n", coord3[0]);
  printf("y vertex 3 %f\n\n", coord3[1]);
  */
  // Generate Jacobian
  double J[2][2];
  J[0][0] = coord1[0] - coord3[0];
  J[0][1] = coord1[1] - coord3[1];
  J[1][0] = coord2[0] - coord3[0];
  J[1][1] = coord2[1] - coord3[1];
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      //printf("J i %d j %d    %f \n",i,j,J[i][j]);
    }
  }
  // Get JJ Matrix
  double JJ[2][2] = {};
  // J inverse
  double Jin[2][2] = {};
  double J1[2][2] = {};
  double J2[2][2] = {};
  Jin[0][0] = J[1][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
  Jin[1][1] = J[0][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
  Jin[1][0] = -J[1][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
  Jin[0][1] = -J[0][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      //printf("Jin %f \n", Jin[i][j]);
    }
  }
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      J1[i][j] = Jin[0][i]*Jin[0][j];
      //printf("J1 %f \n", J1[i][j]);
    }
  }
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      J2[i][j] = Jin[1][i]*Jin[1][j];
      //printf("J2 %f \n", J2[i][j]);
    }
  }
  for (int i = 0; i < 2; i++){
    for (int j = 0; j < 2; j++){
      JJ[i][j] = (J1[i][j]+J2[i][j])*(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      //printf("JJ %f \n", JJ[i][j]);
    }
  }

  // Create Del Matrices and Multiply
  double del[3][2] = {{1,0},{0,1},{-1,-1}};
  double delT[2][3] = {{1,0,-1},{0,1,-1}};
  double delJJ[3][2] = {};
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 2; j++){
      for (int k = 0; k < 2; k++){
        delJJ[i][j] += del[i][k]*JJ[k][j];
        //printf("delJJ i %d j %d    %f \n",i,j,delJJ[i][j]);
      }
    }
  }
  // Create complete matrix, and also do integral (divide by 2)
  double delJJdelT[3][3] = {};
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      for (int k = 0; k < 2; k++){
        delJJdelT[i][j] += delJJ[i][k]*delT[k][j]/2;
      }
    }
  }

  // Generate the forcing vector
  //double f[3] = {};




  // Now simply assemble the contribution and push it back
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      contribution c;
      c.coefficient = delJJdelT[i][j];
      c.known = 0;
      //c.row = pumi_ment_getID(adjacent[i]);
      //c.column = pumi_ment_getID(adjacent[j]);
      c.row = pumi_node_getNumber (numbering, adjacent[i]);
      c.column = pumi_node_getNumber (numbering, adjacent[j]);
      region_contributions.push_back(c);
      //printf("Row %d \n", c.row);
      //printf("Column%d \n", c.column);
      //printf("contribution coefficient %f \n", c.coefficient);
    }
  }
}

// Element routine for T6
void T6(pMeshEnt e, std::vector<contribution> &region_contributions, pNumbering numbering){
  // Get coordiantes of vertices
  int id = pumi_ment_getID(e);
  //printf("\n\n\n\nRegion %d:\n\n",id);
  Adjacent adjacentv;
  pumi_ment_getAdjacent(e,0,adjacentv);
  double coord1[3] = {};
  double coord2[3] = {};
  double coord3[3] = {};
  pumi_node_getCoord(adjacentv[0], 0, coord1);
  pumi_node_getCoord(adjacentv[1], 0, coord2);
  pumi_node_getCoord(adjacentv[2], 0, coord3);
  /*
  printf("x vertex 1 %f\n", coord1[0]);
  printf("y vertex 1 %f\n", coord1[1]);
  printf("x vertex 2 %f\n", coord2[0]);
  printf("y vertex 2 %f\n", coord2[1]);
  printf("x vertex 3 %f\n", coord3[0]);
  printf("y vertex 3 %f\n", coord3[1]);
  */
  // Get coordinates of the edge nodes in proper ordering
  double coord4[3] = {};
  double coord5[3] = {};
  double coord6[3] = {};
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
  /*
  printf("x vertex 4 %f\n", coord4[0]);
  printf("y vertex 4 %f\n", coord4[1]);
  printf("x vertex 5 %f\n", coord5[0]);
  printf("y vertex 5 %f\n", coord5[1]);
  printf("x vertex 6 %f\n", coord6[0]);
  printf("y vertex 6 %f\n\n", coord6[1]);
  */

  // Create x and y vectors
  //double xeT[6] = {coord1[0],coord2[0],coord3[0],coord4[0],coord5[0],coord6[0]};
  //double yeT[6] = {coord1[1],coord2[1],coord3[1],coord4[1],coord5[1],coord6[1]};
  // Create s and t vectors
  // no need for these
  //double s[6] = {0,1,0,0.5,0.5,0};
  //double t[6] = {0,0,1,0,0.5,0.5};
  // Matrix J depends on the quadrature rule;
  // Matrices del and delT depend on the quadrature rules
  // Use 3 point triangle quadrature
  double Sin[3] = {2.0/3,1.0/6,1.0/6};
  double tin[3] = {1.0/6,2.0/3,1.0/6};
  double Win[3] = {1.0/3,1.0/3,1.0/3};

  // create analytical del matrix (split into the three components
  //         constant, s coefficient, and t coefficient)
  double del_analytical[6][2][3] = { {{-3,4,4},{-3,4,4}},
  {{-1,4,0},{0,0,0}},
  {{0,0,0},{-1,0,4}},
  {{4,-8,-4},{0,-4,0}},
  {{0,0,4},{0,4,0}},
  {{0,0,-4},{4,-4,-8}} };
  double result[6][6] = {};
  // calculate the numerical del matrix etc. using the three point quadrature
  // Loop for the three points and then sum up the resulting 6*6 matrix
  for (int l = 0; l < 3; l++){
    //printf("Quadrature %d\n", l);
    double del[6][2] = {};
    for (int i = 0; i < 6; i++){
      for (int j = 0; j < 2; j++){
        // Calculate matrix for this particular quadrature point
        del[i][j] = del_analytical[i][j][0] + del_analytical[i][j][1]*Sin[l] + del_analytical[i][j][2]*tin[l];
        //printf("del matrix i %d j %d    %f\n", i,j,del[i][j]);
        //printf("del analytical coefficients %f %f %f\n",del_analytical[i][j][0],del_analytical[i][j][1],del_analytical[i][j][2]);
        //printf("quadrature points %f %f %f \n", Sin[l], tin[l], Win[l]);
      }
    }

    for (int i = 0; i < 6; i++){
      //printf("del matrix   ");
      for (int j = 0; j < 2; j++){
        //printf(" %f ",del[i][j]);
      }
      //printf("\n");
    }


    // transpose del matrix
    double delT[2][6] = {};
    for (int i = 0; i < 6; i++){
      for (int j = 0; j < 2; j++){
        delT[j][i] = del[i][j];
        //printf("del matrix i %d j %d    %f\n", i,j,del[i][j]);
      }
    }

    /*
    for (int i = 0; i < 2; i++){
    printf("delT matrix   ");
    for (int j = 0; j < 6; j++){
    printf(" %f ",delT[i][j]);
  }
  printf("\n");
}
*/


// Generate Jacobian
// The edges are linear (subparametric) so simply use the T3 Jacobian!
double J[2][2];
J[0][0] = coord1[0] - coord3[0];
J[0][1] = coord1[1] - coord3[1];
J[1][0] = coord2[0] - coord3[0];
J[1][1] = coord2[1] - coord3[1];
for (int i = 0; i < 2; i++){
  for (int j = 0; j < 2; j++){
    //printf("J i %d j %d    %f \n",i,j,J[i][j]);
  }
}
// Get JJ Matrix
double JJ[2][2] = {};
// J inverse
double Jin[2][2] = {};
double J1[2][2] = {};
double J2[2][2] = {};
Jin[0][0] = J[1][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
Jin[1][1] = J[0][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
Jin[1][0] = -J[1][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
Jin[0][1] = -J[0][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
for (int i = 0; i < 2; i++){
  for (int j = 0; j < 2; j++){
    //printf("Jin %f \n", Jin[i][j]);
  }
}
for (int i = 0; i < 2; i++){
  for (int j = 0; j < 2; j++){
    J1[i][j] = Jin[0][i]*Jin[0][j];
    //printf("J1 %f \n", J1[i][j]);
  }
}
for (int i = 0; i < 2; i++){
  for (int j = 0; j < 2; j++){
    J2[i][j] = Jin[1][i]*Jin[1][j];
    //printf("J2 %f \n", J2[i][j]);
  }
}
for (int i = 0; i < 2; i++){
  for (int j = 0; j < 2; j++){
    JJ[i][j] = (J1[i][j]+J2[i][j])*(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
    //printf("JJ %f \n", JJ[i][j]);
  }
}





double delJJ[6][2] = {};
for (int i = 0; i < 6; i++){
  for (int j = 0; j < 2; j++){
    for (int k = 0; k < 2; k++){
      delJJ[i][j] += del[i][k]*JJ[k][j];
      //printf("delJJ i %d j %d    %f \n",i,j,delJJ[i][j]);
    }
  }
}


/*
for (int i = 0; i < 6; i++){
printf("delJJ matrix   ");
for (int j = 0; j < 2; j++){
printf(" %f ",delJJ[i][j]);
}
printf("\n");
}
*/


// Create complete matrix, and also do integral (divide by 2) and use weight
//  (taking advantage of the common weight for all of them!!!!)
// Also multiply by the Jacobian
double delJJdelT[6][6] = {};
for (int i = 0; i < 6; i++){
  for (int j = 0; j < 6; j++){
    for (int k = 0; k < 2; k++){
      delJJdelT[i][j] += delJJ[i][k]*delT[k][j];
    }
  }
}


for (int i = 0; i < 6; i++){
  //printf("delJJdelT matrix   ");
  for (int j = 0; j < 6; j++){
    //printf(" %f ",delJJdelT[i][j]);
  }
  //printf("\n");
}


// Add all this to the result
//printf("Adding to Result\n");
for (int i = 0; i < 6; i++){
  for (int j = 0; j < 6; j++){
    result[i][j] += delJJdelT[i][j];
  }
}
}
for (int i = 0; i < 6; i++){
  //printf("Result matrix   ");
  for (int j = 0; j < 6; j++){
    //printf(" %f ",result[i][j]);
  }
  //printf("\n");
}
// do the weight and jacobian business here
// Now send out the assembly
for (int i = 0; i < 6; i++){
  for (int j = 0; j < 6; j++){
    contribution c;
    c.coefficient = result[i][j];
    c.known = 0;
    // Need to get the proper IDs of the nodes!!!
    if (i<3){
      //c.row = pumi_ment_getID(adjacentv[i]);
      c.row = pumi_node_getNumber (numbering, adjacentv[i]);
    }
    else{
      //c.row = pumi_ment_getID(sharededges[i-3]);
      c.row = pumi_node_getNumber (numbering, sharededges[i-3]);
    }
    if (j<3){
      c.column = pumi_node_getNumber (numbering, adjacentv[j]);
      //pumi_node_getNumber (numbering, adjacentv[j]);
    }
    else{
      //c.column = pumi_ment_getID(sharededges[j-3]);
      c.column = pumi_node_getNumber (numbering, sharededges[j-3]);
    }
    //c.row = pumi_ment_getID(adjacent[i]);
    //c.column = pumi_ment_getID(adjacent[j]);
    region_contributions.push_back(c);
    //printf("Row %d \n", c.row);
    //printf("Column%d \n", c.column);
    //printf("contribution coefficient %f \n", c.coefficient);
    // Verify This !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // fix the determinant division multiplication business
  }
}
}

// Element routine for Q4
void Q4(pMeshEnt e, std::vector<contribution> &region_contributions, pNumbering numbering){
  int id = pumi_ment_getID(e);
  //printf("\n\n\n\nRegion %d:\n\n",id);
  // Get global coordinates
  Adjacent adjacent;
  pumi_ment_getAdjacent(e,0,adjacent);
  double coord1[3];
  double coord2[3];
  double coord3[3];
  double coord4[3];
  pumi_node_getCoord(adjacent[0], 0, coord1);
  pumi_node_getCoord(adjacent[1], 0, coord2);
  pumi_node_getCoord(adjacent[2], 0, coord3);
  pumi_node_getCoord(adjacent[3], 0, coord4);
  /*
  printf("x vertex 1 %f\n", coord1[0]);
  printf("y vertex 1 %f\n", coord1[1]);
  printf("x vertex 2 %f\n", coord2[0]);
  printf("y vertex 2 %f\n", coord2[1]);
  printf("x vertex 3 %f\n", coord3[0]);
  printf("y vertex 3 %f\n", coord3[1]);
  printf("x vertex 4 %f\n", coord4[0]);
  printf("y vertex 4 %f\n\n", coord4[1]);
  */
  // Generate xeT and yeT vectors
  //double xeT[4] = {coord1[0],coord2[0],coord3[0],coord4[0]};
  //double yeT[4] = {coord1[1],coord2[1],coord3[1],coord4[1]};


  // Quadrature Rules
  // Bickford pg. 105
  // Use 2 point quadrature
  double Wi[2] = {1.0,1.0};
  double ui[2] = {0.577350269,-0.577350269};
  double result[4][4] = {};


  // Note: two loops required for quadrature
  // Use s and t as loops for clarity
  for (int s = 0; s < 2; s++){
    for (int t = 0; t < 2; t++){


      // Matrices to generate at the quarature points:
      // Jacobian
      // del


      // pg 299 Bickford
      double J[2][2] = {};
      double J_analytical[2][2][3] = {{{((coord2[0]-coord1[0])+(coord3[0]-coord4[0]))/4,0,(-(coord2[0]-coord1[0])+(coord3[0]-coord4[0]))/4},
      {((coord2[1]-coord1[1])+(coord3[1]-coord4[1]))/4,0,(-(coord2[1]-coord1[1])+(coord3[1]-coord4[1]))/4}},
      {{((coord3[0]-coord2[0])+(coord4[0]-coord1[0]))/4,((coord3[0]-coord2[0])-(coord4[0]-coord1[0]))/4,0},
      {((coord3[1]-coord2[1])+(coord4[1]-coord1[1]))/4,((coord3[1]-coord2[1])-(coord4[1]-coord1[1]))/4,0}}};
      // Bickford pg 309
      double del[4][2] = {};
      double del_analytical[4][2][3] = {{{-1.0/4,0,1.0/4},{-1.0/4,1.0/4,0}},
      {{1.0/4,0,-1.0/4},{-1.0/4,-1.0/4,0}},
      {{1.0/4,0,1.0/4},{1.0/4,1.0/4,0}},
      {{-1.0/4,0,-1.0/4},{1.0/4,-1.0/4,0}}};

      // Calculate J and del at the quadrature point
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          J[i][j] = J_analytical[i][j][0] + ui[s]*J_analytical[i][j][1] + ui[t]*J_analytical[i][j][2];
        }
      }
      for (int i = 0; i < 4; i++){
        for (int j = 0; j < 2; j++){
          del[i][j] = del_analytical[i][j][0] + ui[s]*del_analytical[i][j][1] + ui[t]*del_analytical[i][j][2];
        }
      }



      // FIX ALL INDICES FOR MULTIPLICATION ETC
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          //printf("J %f \n", J[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          //printf("J_analytical 1 %f \n", J_analytical[i][j][0]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          //printf("J_analytical 2 %f \n", J_analytical[i][j][1]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          //printf("J_analytical 3 %f \n", J_analytical[i][j][2]);
        }
      }



      // Get JJ Matrix
      double JJ[2][2] = {};
      double Jin[2][2] = {};
      double J1[2][2] = {};
      double J2[2][2] = {};
      Jin[0][0] = J[1][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      Jin[1][1] = J[0][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      Jin[1][0] = -J[1][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      Jin[0][1] = -J[0][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          //printf("Jin %f \n", Jin[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          J1[i][j] = Jin[0][i]*Jin[0][j];
          //printf("J1 %f \n", J1[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          J2[i][j] = Jin[1][i]*Jin[1][j];
          //printf("J2 %f \n", J2[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          JJ[i][j] = (J1[i][j]+J2[i][j])*(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
          //printf("JJ %f \n", JJ[i][j]);
        }
      }

      // Create delT matrix
      double delT[2][4] = {};
      for (int i = 0; i < 4; i++){
        for (int j = 0; j < 2; j++){
          delT[j][i] = del[i][j];
          //printf("delT matrix i %d j %d    %f\n", i,j,del[i][j]);
        }
      }


      double delJJ[4][2] = {};
      for (int i = 0; i < 4; i++){
        for (int j = 0; j < 2; j++){
          for (int k = 0; k < 2; k++){
            delJJ[i][j] += del[i][k]*JJ[k][j];
            //printf("delJJ i %d j %d    %f \n",i,j,delJJ[i][j]);
          }
        }
      }
      // Create matrix for quadrature point
      double delJJdelT[4][4] = {};
      for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
          for (int k = 0; k < 2; k++){
            delJJdelT[i][j] += delJJ[i][k]*delT[k][j];
          }
        }
      }
      for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
          result[i][j] += delJJdelT[i][j]*Wi[s]*Wi[t];
        }
      }
      // add to result matrix

    }
  }



  // Now simply assemble the contribution and push it back
  for (int i = 0; i < 4; i++){
    for (int j = 0; j < 4; j++){
      contribution c;
      c.coefficient = result[i][j];
      c.known = 0;
      //c.row = pumi_ment_getID(adjacent[i]);
      //c.column = pumi_ment_getID(adjacent[j]);
      c.row = pumi_node_getNumber (numbering, adjacent[i]);
      c.column = pumi_node_getNumber (numbering, adjacent[j]);
      region_contributions.push_back(c);
      /*
      printf("Row %d \n", c.row);
      printf("Column%d \n", c.column);
      printf("contribution coefficient %f \n", c.coefficient);
      */
    }
  }
}

// Element routine for Q8
void Q8(pMeshEnt e, std::vector<contribution> &region_contributions, pNumbering numbering){
  int id = pumi_ment_getID(e);
  //printf("\n\n\n\nRegion %d:\n\n",id);
  // Get global coordinates
  Adjacent adjacentv;
  pumi_ment_getAdjacent(e,0,adjacentv);
  double coord1[3];
  double coord2[3];
  double coord3[3];
  double coord4[3];
  pumi_node_getCoord(adjacentv[0], 0, coord1);
  pumi_node_getCoord(adjacentv[1], 0, coord2);
  pumi_node_getCoord(adjacentv[2], 0, coord3);
  pumi_node_getCoord(adjacentv[3], 0, coord4);
  /*
  printf("x vertex 1 %f\n", coord1[0]);
  printf("y vertex 1 %f\n", coord1[1]);
  printf("x vertex 2 %f\n", coord2[0]);
  printf("y vertex 2 %f\n", coord2[1]);
  printf("x vertex 3 %f\n", coord3[0]);
  printf("y vertex 3 %f\n", coord3[1]);
  printf("x vertex 4 %f\n", coord4[0]);
  printf("y vertex 4 %f\n\n", coord4[1]);
  */
  // Generate xeT and yeT vectors
  //double xeT[4] = {coord1[0],coord2[0],coord3[0],coord4[0]};
  //double yeT[4] = {coord1[1],coord2[1],coord3[1],coord4[1]};




  ////////////////////////////////////////////////////////////////////
  // Get the rest of the vertices using the T6 algorithm
  // Get coordinates of the edge nodes in proper ordering
  double coord5[3] = {};
  double coord6[3] = {};
  double coord7[3] = {};
  double coord8[3] = {};
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
      pumi_node_getCoord(sharededge, 0, coord5);
      break;
      case 1:
      pumi_node_getCoord(sharededge, 0, coord6);
      break;
      case 2:
      pumi_node_getCoord(sharededge, 0, coord7);
      break;
      case 3:
      pumi_node_getCoord(sharededge, 0, coord8);
      break;
      default:
      printf("Unknown Error!!!!\n");
    }


  }/////////////////////////////////////////
  /*
  printf("x vertex 5 %f\n", coord5[0]);
  printf("y vertex 5 %f\n", coord5[1]);
  printf("x vertex 6 %f\n", coord6[0]);
  printf("y vertex 6 %f\n", coord6[1]);
  printf("x vertex 7 %f\n", coord7[0]);
  printf("y vertex 7 %f\n", coord7[1]);
  printf("x vertex 8 %f\n", coord8[0]);
  printf("y vertex 8 %f\n\n", coord8[1]);
  */



  // Quadrature Rules
  // Bickford pg. 105
  // Use 2 point quadrature
  double Wi[2] = {1.0,1.0};
  double ui[2] = {0.577350269,-0.577350269};
  double result[8][8] = {};


  // Note: two loops required for quadrature
  // Use s and t as loops for clarity
  for (int s = 0; s < 2; s++){
    for (int t = 0; t < 2; t++){


      // Matrices to generate at the quarature points:
      // Jacobian
      // del

      // Jacobian is the same as the Q4 element since the elements are subparametric and have linear edges


      // pg 299 Bickford
      double J[2][2] = {};
      double J_analytical[2][2][3] = {{{((coord2[0]-coord1[0])+(coord3[0]-coord4[0]))/4,0,(-(coord2[0]-coord1[0])+(coord3[0]-coord4[0]))/4},
      {((coord2[1]-coord1[1])+(coord3[1]-coord4[1]))/4,0,(-(coord2[1]-coord1[1])+(coord3[1]-coord4[1]))/4}},
      {{((coord3[0]-coord2[0])+(coord4[0]-coord1[0]))/4,((coord3[0]-coord2[0])-(coord4[0]-coord1[0]))/4,0},
      {((coord3[1]-coord2[1])+(coord4[1]-coord1[1]))/4,((coord3[1]-coord2[1])-(coord4[1]-coord1[1]))/4,0}}};
      // Bickford pg 309
      double del[8][2] = {  {(1-ui[t])*(2*ui[s]+ui[t])/4 , (1-ui[s])*(ui[s]+2*ui[t])/4},
      {(1-ui[t])*(2*ui[s]-ui[t])/4 , (1+ui[s])*(-ui[s]+2*ui[t])/4},
      {(1+ui[t])*(2*ui[s]+ui[t])/4 , (1+ui[s])*(ui[s]+2*ui[t])/4},
      {(1+ui[t])*(2*ui[s]-ui[t])/4 , (1-ui[s])*(-ui[s]+2*ui[t])/4},
      {-ui[s]*(1-ui[t]) , -(1-ui[s]*ui[s])/2},
      {(1-ui[t]*ui[t])/2 , -ui[t]*(1+ui[s])},
      {-ui[s]*(1+ui[t]) , (1-ui[s]*ui[s])/2},
      {-(1-ui[t]*ui[t])/2 , -ui[t]*(1-ui[s])}  };

      // Calculate J and del at the quadrature point
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          J[i][j] = J_analytical[i][j][0] + ui[s]*J_analytical[i][j][1] + ui[t]*J_analytical[i][j][2];
        }
      }






      // Get JJ Matrix
      double JJ[2][2] = {};
      double Jin[2][2] = {};
      double J1[2][2] = {};
      double J2[2][2] = {};
      Jin[0][0] = J[1][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      Jin[1][1] = J[0][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      Jin[1][0] = -J[1][0]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      Jin[0][1] = -J[0][1]/(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          //printf("Jin %f \n", Jin[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          J1[i][j] = Jin[0][i]*Jin[0][j];
          //printf("J1 %f \n", J1[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          J2[i][j] = Jin[1][i]*Jin[1][j];
          //printf("J2 %f \n", J2[i][j]);
        }
      }
      for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
          JJ[i][j] = (J1[i][j]+J2[i][j])*(J[1][1]*J[0][0]-J[1][0]*J[0][1]);
          //printf("JJ %f \n", JJ[i][j]);
        }
      }

      // Create delT matrix
      double delT[2][8] = {};
      for (int i = 0; i < 8; i++){
        for (int j = 0; j < 2; j++){
          delT[j][i] = del[i][j];
          //printf("delT matrix i %d j %d    %f\n", i,j,del[i][j]);
        }
      }


      double delJJ[8][2] = {};
      for (int i = 0; i < 8; i++){
        for (int j = 0; j < 2; j++){
          for (int k = 0; k < 2; k++){
            delJJ[i][j] += del[i][k]*JJ[k][j];
            //printf("delJJ i %d j %d    %f \n",i,j,delJJ[i][j]);
          }
        }
      }
      // Create matrix for quadrature point
      double delJJdelT[8][8] = {};
      for (int i = 0; i < 8; i++){
        for (int j = 0; j < 8; j++){
          for (int k = 0; k < 2; k++){
            delJJdelT[i][j] += delJJ[i][k]*delT[k][j];
          }
        }
      }
      for (int i = 0; i < 8; i++){
        for (int j = 0; j < 8; j++){
          result[i][j] += delJJdelT[i][j]*Wi[s]*Wi[t];
        }
      }
      // add to result matrix
      // print result matrix
      for (int i = 0; i < 8; i++){
        //printf("result     %d", i);
        for (int j = 0; j < 8; j++){
          //printf(" %2.2f ", result[i][j]);
        }
        //printf("\n");
      }
    }
  }



  // Now simply assemble the contribution and push it back
  for (int i = 0; i < 8; i++){
    for (int j = 0; j < 8; j++){
      contribution c;
      c.coefficient = result[i][j];
      c.known = 0;
      // Need to get the proper IDs of the nodes!!!
      if (i<4){
        //c.row = pumi_ment_getID(adjacentv[i]);
        c.row = pumi_node_getNumber (numbering, adjacentv[i]);
      }
      else{
        //c.row = pumi_ment_getID(sharededges[i-3]);
        c.row = pumi_node_getNumber (numbering, sharededges[i-4]);
      }
      if (j<4){
        c.column = pumi_node_getNumber (numbering, adjacentv[j]);
        //pumi_node_getNumber (numbering, adjacentv[j]);
      }
      else{
        //c.column = pumi_ment_getID(sharededges[j-3]);
        c.column = pumi_node_getNumber (numbering, sharededges[j-4]);
      }
      //c.row = pumi_ment_getID(adjacent[i]);
      //c.column = pumi_ment_getID(adjacent[j]);
      region_contributions.push_back(c);
      /*
      printf("Row %d \n", c.row);
      printf("Column%d \n", c.column);
      printf("contribution coefficient %f \n", c.coefficient);
      */
    }
  }
}

// Generic region routine, will region contributions for all defined elements
contribution region_routine(pMesh mesh, pMeshEnt e, pNumbering numbering, std::vector<contribution> &region_contributions){
  //contribution region_contribution;
  switch (get_element_type(mesh,e)) {
    case 0:
    printf("      Unrecognised element!!!!!!!!!!!\n");
    break;
    case 1:
    //printf("T3\n");
    T3(e, region_contributions, numbering);
    break;
    case 2:
    //printf("T6\n");
    T6(e, region_contributions, numbering);
    break;
    case 3:
    //printf("Q4\n");
    Q4(e, region_contributions, numbering);
    break;
    case 4:
    //printf("Q8\n");
    Q8(e, region_contributions, numbering);
    break;
    default:
    printf("      Unknown Error in Element routine!!!!!!!!!!\n");
    break;
  }
}

// Also implement the forcing contribution


// Now need to implement the edge boundary conditions and the essential boundary conditon enforcer
// This will be done using the boundary conditions vectors that we have created.


// Edge routine for boundaries
void edge_routine(pMesh mesh, boundary_struct edge, pNumbering numbering, std::vector<contribution> &edge_contributions){
  //Hughes pg. 109
  // Get the boundary number using the already defined function
  // Apply the proper quadrature at that boundary and create the contributions
  // Boundary may have nodes on it too, so implement the seperate quadrature for it also.
  // Detect the proper type of the edge and simply use the proper procedure accordingly

  // Get the length of the edge and then transform it
  // Get the coordinates of the adjacent vertices.
  Adjacent adjacent;
  pumi_ment_getAdjacent(edge.e,0,adjacent);
  double coord1[3];
  double coord2[3];
  pumi_node_getCoord(adjacent[0], 0, coord1);
  pumi_node_getCoord(adjacent[1], 0, coord2);
  // Calculate the length
  double len = sqrt( pow(coord1[0]-coord2[0],2) + pow(coord1[1]-coord2[1],2) );



  if (!hasNode(mesh,edge.e)){
    // Use linear function quadrature to get the edge contribution
    // Bickford pg.256-7
    double alpha = edge.first;
    double a[2][2] = {};
    a[0][0] = len/12*(alpha*3 + alpha);
    a[0][1] = len/12*(alpha + alpha);
    a[1][0] = len/12*(alpha + alpha);
    a[1][1] = len/12*(alpha + 3*alpha);

    double he = edge.second;
    double h[2] = {};
    h[0] = len/6*(2*he+he);
    h[1] = len/6*(he+2*he);
    // Push back right here to remain in scope
    for (int i = 0; i < 2; i++){
      for (int j = 0; j < 2; j++){
        contribution c;
        c.coefficient = a[i][j];
        c.known = 0;
        //c.row = pumi_ment_getID(adjacent[i]);
        //c.column = pumi_ment_getID(adjacent[j]);
        c.row = pumi_node_getNumber (numbering, adjacent[i]);
        c.column = pumi_node_getNumber (numbering, adjacent[j]);
        edge_contributions.push_back(c);
        //printf("Row %d ", c.row);
        //printf("Column %d ", c.column);
        //printf("contribution coefficient %f \n", c.coefficient);
      }
    }
    // push back coefficients
    for (int i = 0; i < 2; i++){
      contribution c;
      c.coefficient = 0;
      c.known = h[i];
      //c.row = pumi_ment_getID(adjacent[i]);
      //c.column = pumi_ment_getID(adjacent[j]);
      c.row = pumi_node_getNumber (numbering, adjacent[i]);
      c.column = 0;
      edge_contributions.push_back(c);
      //printf("Row %d ", c.row);
      //printf("known %f \n", c.known);
    }



  }
  else {
    // Use quadratic function quadrature
    // Bickford pg 331
    double alpha = edge.first;
    double a[3][3] = {};
    a[0][0] = len/420*(alpha*39 + alpha*20 - alpha*3);
    a[0][1] = len/420*(alpha*20 + alpha*16 - alpha*8);
    a[0][2] = len/420*(-alpha*3 -alpha*8 - alpha*3);

    a[1][0] = len/420*(alpha*20 + alpha*16 - alpha*8);
    a[1][1] = len/420*(alpha*16 + alpha*192 + alpha*16);
    a[1][2] = len/420*(-alpha*8 + alpha*16 + alpha*20);

    a[2][0] = len/420*(-alpha*3 -alpha*8 - alpha*3);
    a[2][1] = len/420*(-alpha*8 + alpha*16 + alpha*20);
    a[2][2] = len/420*(-alpha*3 + alpha*20 + alpha*39);

    double he = edge.second;
    double h[3] = {};
    h[0] = len/30*(4*he+2*he-he);
    h[1] = len/30*(2*he+16*he+2*he);
    h[2] = len/30*(-he+2*he+4*he);

    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        contribution c;
        c.coefficient = a[i][j];
        c.known = 0;
        //c.row = pumi_ment_getID(adjacent[i]);
        //c.column = pumi_ment_getID(adjacent[j]);
        if (1==i){
          c.row = pumi_node_getNumber (numbering, edge.e);
        }
        else if (2==i){
          c.row = pumi_node_getNumber (numbering, adjacent[i-1]);
        }
        else {
          c.row = pumi_node_getNumber (numbering, adjacent[i]);
        }
        if (1==j){
          c.column = pumi_node_getNumber (numbering, edge.e);
        }
        else if(2==j){
          c.column = pumi_node_getNumber (numbering, adjacent[j-1]);
        }
        else {
          c.column = pumi_node_getNumber (numbering, adjacent[j]);
        }
        //c.column = pumi_node_getNumber (numbering, adjacent[j]);
        edge_contributions.push_back(c);
        //printf("Row %d ", c.row);
        //printf("Column %d ", c.column);
        //printf("contribution coefficient %f \n", c.coefficient);
      }
    }

    for (int i = 0; i < 3; i++){
      contribution c;
      c.coefficient = 0;
      c.known = h[i];
      //c.row = pumi_ment_getID(adjacent[i]);
      //c.column = pumi_ment_getID(adjacent[j]);
      if (1==i){
        c.row = pumi_node_getNumber (numbering, edge.e);
      }
      else if (2==i){
        c.row = pumi_node_getNumber (numbering, adjacent[i-1]);
      }
      else {
        c.row = pumi_node_getNumber (numbering, adjacent[i]);
      }
      c.column = 0;

      //c.column = pumi_node_getNumber (numbering, adjacent[j]);
      edge_contributions.push_back(c);
      //printf("Row %d ", c.row);
      //printf("known %f \n", c.known);
      //printf("Row %d \n", c.row);
      //printf("Column%d \n", c.column);
      //printf("contribution coefficient %f \n", c.coefficient);
    }
  }
}



// Force Contributor (face)

// Constrain a node



//
