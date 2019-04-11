#include <PCU.h>
#include <pumi.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cmath>
using namespace std;
#include "reorder.cc"
#include "element_routines.cc"
#include "matrix_assembly.cc"
#include "postprocessing_routines.cc"
//#include <petscksp.h>


//double boundary_conditions(pMeshEnt e);

//std::pair<bool,double> essential_BC(int boundary_number);
//std::pair<double,double> natural_BC(int boundary_number);

int main(int argc, char** argv)
{

  //////////////////// Startup code. ////////////////

  // Check if files have been passed properly
  if (argc != 3) {
    printf("usage: %s dmg_file.dmg smb_file_without_0.smb\n", argv[0]);
    return 0;
  }

  // Load the files passed
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load(argv[1],"mesh");
  pMesh mesh = pumi_mesh_load(geom,argv[2],1);

  // Convert some mesh to Lagrange or Serendipity
  if(!strcmp (argv[1], "reorder_a.dmg")) {
    pumi_mesh_setShape(mesh,pumi_shape_getSerendipity());
  }
  else{
    pumi_mesh_setShape(mesh,pumi_shape_getLagrange(2));
  }
  pumi_mesh_print(mesh);



  // Inform that files have successfully been loaded
  printf("\nMesh and Geometry Successfully loaded\n");
  // Write this mesh file for debugging purposes
  pumi_mesh_write(mesh,"originalmesh","vtk");


  // Now that the mesh files have successfully been loaded, reorder them
  pNumbering numbering;
  numbering = reorder_mesh(mesh,geom);
  cout << "Completed reordering\n";
  pumi_mesh_write(mesh,argv[2],"vtk");
  cout << "Mesh file written to file\n";






  // Check the element type routine, if it is working
  pMeshEnt e;
  pMeshIter it = mesh->begin(2);
  while ((e = mesh->iterate(it))){
    int element_type = get_element_type(mesh,e);
    //printf("The element type is %d\n", element_type);
  }
  mesh->end(it);






  // Get all the boundary edges
  std::vector<boundary_struct> boundary_edges;
  get_all_boundary_edges(geom,mesh, boundary_edges);
  // Print the boundary numbers
  for (std::vector<boundary_struct>::iterator it = boundary_edges.begin(); it!= boundary_edges.end(); ++it){
    boundary_struct this_edge = *it;
    //printf("Edge is on face %d \n", this_edge.boundary);
  }

  // Get all the boundary vertices
  std::vector<boundary_struct> boundary_verts;
  get_all_boundary_nodes(mesh, boundary_edges, boundary_verts, numbering);
  // Now print all the boundary nodes
  for (std::vector<boundary_struct>::iterator it = boundary_verts.begin(); it!= boundary_verts.end(); ++it){
    boundary_struct this_vert = *it;
    //printf("Vertex %d is on %d \n", pumi_node_getNumber(numbering, this_vert.e), this_vert.boundary);
  }
  // Try to find which boundary something is on
  std::vector<int> list;
  int tofind = 0;
  pMeshEnt ment = pumi_mesh_findEnt(mesh,0,tofind);
  get_bound_num(ment, boundary_verts, list);
  for (int i = 0; i< list.size(); i++){
    //printf("The vertex %d is on %d\n",tofind,list[i]);
  }


  // Check is redordering is required
  it = mesh->begin(2);
  while ((e = mesh->iterate(it))){
    Adjacent adjacent;
    pumi_ment_getAdjacent(e,0,adjacent);
    if(reorder_verts(adjacent)){
      printf("                   Vertices need reordering \n");
      return 0;
    }
  }
  mesh->end(it);
  printf("Vertices do not need reordering\n");


  // Now we are ready generate the contributions
  std::vector<contribution> all_contributions;
  // Generate the region contributions
  it = mesh->begin(2);
  while ((e = mesh->iterate(it))){
    //Q8(e, all_contributions, numbering);
    region_routine(mesh, e, numbering, all_contributions);
  }
  mesh->end(it);
  printf("Generated region contributions\n");

  // Loop over all the boundary edges to get the contributions

  for (int i = 0; i < boundary_edges.size(); i++){
    edge_routine(mesh, boundary_edges[i], numbering, all_contributions);
  }
  printf("Generated edge contributions \n");




  // Generate the global stiffness matrix and the b vector
  // Generate the vector of knowns
  int m = pumi_numbering_getNumNode(numbering);
  printf("Total number of nodes is %d \n", m);
  double A[m][m] = {};
  double b[m] = {};
  for (int i = 0; i < all_contributions.size(); i++){
    A[all_contributions[i].row][all_contributions[i].column] += all_contributions[i].coefficient;
    //printf("the coefficient is %f \n", all_contributions[i].coefficient);
    b[all_contributions[i].row] += all_contributions[i].known;
  }

  // print the global stiffness matrix
  /*
  for (int i = 0; i < m; i++){
  printf("Row %d ", i);
  for (int j = 0; j < m; j++){
  printf(" %3.2f ", A[i][j]);
}
printf("\n");
}
*/

// Check the symmetry of the system
for (int i = 0; i < m; i++){
  for (int j = 0; j < m; j++){
    //printf(" %.10f \n", A[i][j]-A[j][i]);
    if (A[i][j]-A[j][i] > 0.00000001){
      printf("                              Global Stiffnesss matrix is NOT symmetrical \n");
      printf("row %d column %d \n", i,j);
      return 0;
    }
  }
  //printf("\n");
}


// Check the known vector
for (int i = 0; i < m; i++){
  //printf("row %d   known  %f \n", i, b[i]);
}

// Enforce the dirichlet boundary condition to make node 0 = 1;
A[0][0] = 1;
b[0] = 1;
for (int i = 1; i < m; i++){
  A[0][i] = 0;
}

// print the global stiffness matrix
/*
for (int i = 0; i < m; i++){
printf("Row %d ", i);
for (int j = 0; j < m; j++){
printf(" %3.2f ", A[i][j]);
}
printf("\n");
}
*/

// Check the known vector
for (int i = 0; i < m; i++){
  printf("row %d   known  %f \n", i, b[i]);
}









// Use some library to solve the assembled matrix


// Now run some code to generate the secondary variables and assign them to some field


// Visualize the solution with Paraview directly if possible




pumi_finalize();
MPI_Finalize();
cout << "                                   End of Program\n";
}

/*
// Define the boundary conditions of the problem here
// pass a node or and edhe to get its boundary condition
bool boundary_condition(pMeshEnt e,std::vector<boundary_struct> &boundary_edges,std::vector<boundary_struct> &boundary_verts){
// Check the type of the entity
// If the entity is a vertex iterate, find its boundaries.
// Iterate over the list and perform switch on all the elements, fill out a vector.
// In case the length of the vetor is 2, give an error that two conditions have been defined on the boundary.
// If not, simply fill out the boundary condition
//
// In case the entity is an edge, simply sind its boundary,
// perform a switch, and return the boundary condition of the entity.

// Find the borders that the entity is on
// Fill out the appropriate lists
//  Run the appropriate boundary code for either

if(list.size() != 0){
// If the entity is on a boundary, use the boundary number to perform a switch
// selction and return the boundary condition
}
else{
// If nothing is found, print out the error message and abort
printf("\n ELEMENT IS NOT ON A BOUNDARY!!! \n");
}
}
*/






// Simple definitions of boundaries.

// Define the essential boundary condition using either the local ID or the boundary number
BC essential_BC(int boundary_number, pMeshEnt e, pNumbering numbering){
  // simply return the known
  BC BCe;
  //int localid = pumi_ment_getID(e);
  // get id of the field node using the numbering
  int number = pumi_node_getNumber (numbering, e);

  BCe.first = 1;
  switch (boundary_number) {
    case 0:
    BCe.first = 0;
    BCe.second = 0;
    return BCe;
    default:

    BCe.first = 0;
    BCe.second = 0;
    return BCe;
  }

  switch (number) {
    case 0:
    BCe.first = 0;
    BCe.second = 0;
    return BCe;
    default:

    BCe.first = 0;
    BCe.second = 0;
    return BCe;
  }
}


// Define natural boundary condition using the boundary number
BC natural_BC(int boundary_number){
  // push back alpha and then h (coefficeient and then known)
  BC BCn;
  switch (boundary_number) {
    case 0:
    BCn.first = 1;
    BCn.second = 1;
    return BCn;

    case 6:
    BCn.first = 1;
    BCn.second = 1;
    return BCn;

    default:
    BCn.first = 0;
    BCn.second = 0;
    return BCn;
  }
}



// Old BC Structure
/*
std::pair<bool,double> essential_BC(int boundary_number){
// simply return the known
std::pair<bool,double> BC;
BC.first = true;
switch (boundary_number) {
case 0:
BC.first = false;
BC.second = 0;
return BC;
default:

BC.first = false;
BC.second = 0;
return BC;
}
}

std::pair<double,double> natural_BC(int boundary_number){
// push back alpha and then h (coefficeient and then known)
std::pair<double,double> BC;
switch (boundary_number) {
case 0:
BC.first = 0.5;
BC.second = 0.8;
return BC;

case 1:
BC.first = 0.5;
BC.second = 0.8;
return BC;

default:
BC.first = 0;
BC.second = 0;
return BC;
}
}
*/
