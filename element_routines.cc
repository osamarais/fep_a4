// This file contains all the functions for the element routines

// Header
std::pair<bool,double> essential_BC(int boundary_number);
std::pair<double,double> natural_BC(int boundary_number);

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
    //printf("Unidentified Element!!!\n");
  }
  //printf("The element_type type identified is %d \n", element_type);
  return element_type;
}








// Get all the mesh edges classified on ALL geometric boundaries

struct boundary_struct{
  int boundary;
  pMeshEnt e;
};

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



// The contribution calculator for a region
// Needs:
// Outputs: a vecor of contributions

// Algorithm:
// 1) get the type of the element in order to determine the shape function to be implemented
// 2) get the adjacent vertices edges in order to get the proper nodes to contribute to
// 3) Start getting the contributions
// 4) Do a mapping, calculate the Jacobian.
// 




// The contribution calculator for Edges



// The essential BC enforcer




























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
  // Needs: boundary lists (verties and edges), mesh entity (node), mesh
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
