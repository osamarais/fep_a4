// This file contains all the functions for the element routines

// Minimum requirements are 3 sided and 4 sided Linear and Quaratic Elements

// The type of elements in the mesh must be identified. Functions for that are placed here
/*

Element enumeration is as follows:
1: Linear Tet
2: Quadratic (Serendipity) Tet
3: Linear Quad
4: Quadratic (Serendipity) Quad

*/

int get_element_type(pMesh mesh){
  int element_type = 0;
  int edgenodes = 0;
  int numberofedges = 0;
  // get nummber of nodes on the edges
  pShape myshape = pumi_mesh_getShape(mesh);
  edgenodes = pumi_shape_getNumNode(myshape, 1);
  // get number of nodes on the faces
  // implement this later if needed
  // get any random mesh entity and check the nmber of adjacent edges to get the shape of the element
  pMeshEnt e;
  pMeshIter it = mesh->begin(2);
  while ((e = mesh->iterate(it))){
    numberofedges = pumi_ment_getNumAdj(e,1);
    break;
  }
  mesh->end(it);
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
// pass std::vector<boundary_edge_struct> to the function to get it filled out
struct boundary_edge_struct{
  int boundary;
  pMeshEnt e;
};
struct boundary_vert_struct{
  int boundary;
  pMeshEnt e;
};

void get_all_boundary_edges(pGeom &g,pMesh &mesh, std::vector<boundary_edge_struct> &mesh_ents){
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
      boundary_edge_struct this_edge;
      this_edge.boundary = pumi_gent_getID(*it1);
      this_edge.e = *it2;
      mesh_ents.push_back(this_edge);
    }
  }
}


// get a list of all the vertices on boundaries
void get_all_boundary_nodes(pMesh &mesh, std::vector<boundary_vert_struct> &mesh_ents){
  // Get the adjacent vertices of all the edges.
  // avoid duplication by checking the list and tallying against e and boundary both
  // allow vertex to be classifoed on more than one boundary
  // Iterate over the edges
  for (std::vector<boundary_edge_struct>::iterator it1 = boundary_edges.begin(); it1!= boundary_edges.end(); ++it1){
    boundary_edge_struct this_edge = *it1;
    // Get the adjacent edges
    // go through the container, check if both boundary and e have not been satisfied. if not push back
    Adjacent adjacent;
    pumi_ment_getAdjacent(this_edge.e,0,adjacent);
    for(int i = 0; i< size(adjacent); i ++){
      for (std::vector<boundary_vert_struct>::iterator it2 = mesh_ents.begin(); it2!= mesh_ents.end(); ++it2){
        printf("Finding something here\n");
      }
    }
  }
}




// Function to check which geometric face (boundary) an element is on (zero is none)
int get_edge_bound_num(pMeshEnt ment, std::vector<std::vector<pMeshEnt>> &mesh_ents){
  // do a nested iteration to check if the entity is on a boundary or not
  int boundary_counter = 0;
  for(std::vector<std::vector<pMeshEnt>>::iterator it1 = mesh_ents.begin() ; it1!=mesh_ents.end(); ++it1){
    boundary_counter++;
    std::vector<pMeshEnt> this_group = *it1;
    for(std::vector<pMeshEnt>::iterator it2 = this_group.begin() ; it2!=this_group.end(); ++it2){
      if(*it2 == ment){
        return boundary_counter;
      }
    }
  }
  return 0;
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


























//











/* backup
void get_all_boundary_edges(pGeom &g,pMesh &mesh, std::vector<std::vector<pMeshEnt>> &mesh_ents){
  //printf("getting boundary elements \n");
  // Iterate over all the geometric edges
  // put mesh entities classified on those edges into the container
  for (pGeomIter it = g->begin(1); it!=g->end(1);++it){
    //printf("looping in Geometry\n");
    pGeomEnt ge = *it;
    std::vector<pMeshEnt> entities_on_this;
    // get all the reverse classified entites on this
    pumi_gent_getRevClas(ge, entities_on_this);
    //printf("got the reverse classification\n");
    // push the subcontainer into the main container
    mesh_ents.push_back(entities_on_this);
    //printf("pushed back the subcontainer\n");
  }
}
*/
