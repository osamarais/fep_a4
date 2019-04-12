#include <PCU.h>
#include <pumi.h>
#include <apf.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cmath>
using namespace std;
#include "reorder.cc"
#include "element_routines.cc"
#include "matrix_assembly.cc"
#include "postprocessing_routines.cc"
#include <Eigen/Dense>
using namespace Eigen;

int main(int argc, char** argv)
{
  bool isproblemquadratic = false;
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
  if (isproblemquadratic){
    // Convert some mesh to Lagrange or Serendipity
    if(!strcmp (argv[1], "reorder_a.dmg")) {
      pumi_mesh_setShape(mesh,pumi_shape_getSerendipity());
    }
    else{
      pumi_mesh_setShape(mesh,pumi_shape_getLagrange(2));
    }
  }
  printf("\nMesh and Geometry Successfully loaded\n");
  // Reorder mesh
  pNumbering numbering;
  numbering = reorder_mesh(mesh,geom);
  printf("Completed reordering\n");
  // Get all the boundary edges
  std::vector<boundary_struct> boundary_edges;
  get_all_boundary_edges(geom,mesh, boundary_edges);
  // Print all the faces the boundary edges are on
  for (std::vector<boundary_struct>::iterator it = boundary_edges.begin(); it!= boundary_edges.end(); ++it){
    boundary_struct this_edge = *it;
    //printf("Edge is on face %d \n", this_edge.boundary);
  }
  // Get all the boundary vertices
  std::vector<boundary_struct> boundary_verts;
  get_all_boundary_nodes(mesh, boundary_edges, boundary_verts, numbering);
  // Print all the faces the boundary vertices are on
  for (std::vector<boundary_struct>::iterator it = boundary_verts.begin(); it!= boundary_verts.end(); ++it){
    boundary_struct this_vert = *it;
    printf("Vertex %d is on face %d \n", pumi_node_getNumber(numbering, this_vert.e), this_vert.boundary);
  }
  // Try to find which boundary an entity is on (pass the appropriate list)
  std::vector<int> list;
  int tofind = 0;
  pMeshEnt ment = pumi_mesh_findEnt(mesh,0,tofind);
  get_bound_num(ment, boundary_verts, list);
  for (int i = 0; i< list.size(); i++){
    //printf("The vertex %d is on %d\n",tofind,list[i]);
  }

  // Check if all the vertices are being returned in clockwise order
  // Return error if violated
  pMeshIter it;
  pMeshEnt e;
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


  // Generate all the elemental stiffness matrices
  std::vector<contribution> all_contributions;
  // Generate the region contributions
  it = mesh->begin(2);
  while ((e = mesh->iterate(it))){
    region_routine(mesh, e, numbering, all_contributions);
  }
  mesh->end(it);
  printf("Generated region contributions\n");
  // Generate boundary contributions
  for (int i = 0; i < boundary_edges.size(); i++){
    edge_routine(mesh, boundary_edges[i], numbering, all_contributions);
  }
  printf("Generated edge contributions \n");




  // Generate the global stiffness matrix and the b vector
  // Generate the vector of knowns
  int m = pumi_numbering_getNumNode(numbering);
  printf("Total number of nodes is %d \n", m);
  MatrixXf A = MatrixXf::Zero(m, m);
  VectorXf b = VectorXf::Zero(m);
  VectorXf u = VectorXf::Zero(m);
  // Assemble the Global stiffness matrix and the vector of knowns
  for (int i = 0; i < all_contributions.size(); i++){
    A(all_contributions[i].row,all_contributions[i].column) += all_contributions[i].coefficient;
    //printf("the coefficient is %f \n", all_contributions[i].coefficient);
    b(all_contributions[i].row) += all_contributions[i].known;
  }
  // Check the symmetry of the system
  for (int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      //printf(" %.10f \n", A[i][j]-A[j][i]);
      if (A(i,j)-A(j,i) > 0.00000001){
        printf("                              Global Stiffnesss matrix is NOT symmetrical \n");
        printf("row %d column %d \n", i,j);
        return 0;
      }
      A(i,j) = A(j,i);
    }
    //printf("\n");
  }
  // Check the known vector
  for (int i = 0; i < m; i++){
    //printf("row %d   known  %f \n", i, b[i]);
  }


  // If a mixed boundary condition has not been specified,
  // essential boundary conditions may need to be enforced





  // Enforce the dirichlet boundary conditions

  for (std::vector<boundary_struct>::iterator it = boundary_verts.begin(); it!= boundary_verts.end(); ++it){
    boundary_struct this_vert = *it;
    //printf("Vertex %d is on face %d \n", pumi_node_getNumber(numbering, this_vert.e), this_vert.boundary);
    BC essential = essential_BC(this_edge.boundary, this_edge.e, numbering);
    if(essential.first){
      int rcn = pumi_numbering_getNumNode(numbering);
      A(rcn,rcn) = 1000000000000000;
      b(rcn) = 1000000000000000*BC.second;
    }

  }








  // Solve the system
  u = A.partialPivLu().solve(b);

  // Apply the solution to the field nodes
  pField myfield =  pumi_field_create(mesh,"Primary Solution",1);
  // Apply to vertex nodes
  it = mesh->begin(0);
  while ((e = mesh->iterate(it))){
    double myvar[1] = {u(pumi_node_getNumber(numbering,e))};
    apf::setScalar(myfield, e, 0, myvar[0]);
  }
  mesh->end(it);
  // Apply to edge nodes
  it = mesh->begin(1);
  while ((e = mesh->iterate(it))){
    if (hasNode(mesh, e)){
      double myvar[1] = {u(pumi_node_getNumber(numbering,e))};
      apf::setScalar(myfield, e, 0, myvar[0]);
    }
  }
  mesh->end(it);




  // Freeze field and write out mesh
  pumi_field_freeze (myfield);
  pumi_mesh_write(mesh,argv[2],"vtk");
  printf("Mesh file written to file\n");
  pumi_finalize();
  MPI_Finalize();
  printf("End of Program\n");
}




// Simple definitions of boundaries using faces.

// Define the essential boundary condition using the face it is classified on
BC essential_BC(int boundary_number, pMeshEnt e, pNumbering numbering){
  BC BCe;
  BCe.first = 1;

  // Use this to enforce essential boundary conditions on a face
  switch (boundary_number) {
    case 0:
    BCe.second = 1;
    return BCe;
    case 2:
    BCe.second = -1;
    return BCe;
    case 6:
    BCe.second = 1;
    return BCe;
    case 16:
    BCe.second = -1;
    return BCe;


    default:
    BCe.first = 0;
    BCe.second = 0;
    return BCe;
  }

  // Use this to enforce essential boundary conditions on specific nodes
  /*
  int number = pumi_node_getNumber (numbering, e);
  switch (number) {
    case 0:
    BCe.first = 0;
    BCe.second = 0;
    return BCe;

    default:
    BCe.first = 0;
    BCe.second = 0;
    return BCe;
  */
  }
}


// Define natural boundary condition using the face it is classified on
BC natural_BC(int boundary_number){
  // push back alpha and then h (coefficeient and then known)
  BC BCn;
  switch (boundary_number) {
    /*
    case 0:
    BCn.first = 10;
    BCn.second = 100;
    return BCn;

    case 2:
    BCn.first = 10;
    BCn.second = -100;
    return BCn;


    case 6:
    BCn.first = 10;
    BCn.second = 100;
    return BCn;

    case 16:
    BCn.first = 10;
    BCn.second = -100;
    return BCn;
    */


    default:
    BCn.first = 0;
    BCn.second = 0;
    return BCn;
  }
}
