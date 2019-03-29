#include <PCU.h>
#include <pumi.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
using namespace std;
#include "reorder.cc"
#include "element_routines.cc"
#include "matrix_assembly.cc"
#include "postprocessing_routines.cc"


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

  // Convert some mesh to Lagrange
  if(!strcmp (argv[1], "reorder_c.dmg")) {
    pumi_mesh_setShape(mesh,pumi_shape_getLagrange(2));
    pumi_mesh_print(mesh);
  }



  // Inform that files have successfully been loaded
  printf("\nMesh and Geometry Successfully loaded\n");
  // Write this mesh file for debugging purposes
  pumi_mesh_write(mesh,"originalmesh","vtk");


  // Now that the mesh files have successfully been loaded, reorder them if need be
  char needtoreorder;
  cout << "\nDo you want to reorder the mesh?\nEnter 'y' if so\n      Enter anything else otherwise\n";
  //cin >> needtoreorder;
  //if ('y'==needtoreorder) {
  if (true) {
    cout << "Reordering the mesh\n";
    // Reorder the mesh
    reorder_mesh(mesh,geom);
    cout << "Completed reordering\n";
    pumi_mesh_write(mesh,argv[2],"vtk");
    cout << "Mesh file written to file\n";
  }
  else {
    cout << "Skipped reordering, proceeding to solve\n";
  }






  // Now the element routines must be run over the mesh elements in a loop while assembling the matrix.
  int element_type = get_element_type(mesh);
  printf("The element type is %d\n", element_type);




  // Get all the boundary edges
  std::vector<boundary_struct> boundary_edges;
  get_all_boundary_edges(geom,mesh, boundary_edges);
  // Print the boundary numbers
  for (std::vector<boundary_struct>::iterator it = boundary_edges.begin(); it!= boundary_edges.end(); ++it){
    boundary_struct this_edge = *it;
    printf("Edge is on face %d \n", this_edge.boundary);
  }
  // Get boundary number for





  std::vector<boundary_struct> boundary_verts;
  // Get all the boundary vertices
  get_all_boundary_nodes(mesh, boundary_edges, boundary_verts);
  // Now print all the boundary nodes
  for (std::vector<boundary_struct>::iterator it = boundary_verts.begin(); it!= boundary_verts.end(); ++it){
    boundary_struct this_vert = *it;
    printf("Vertex %d is on %d \n", pumi_ment_getID(this_vert.e), this_vert.boundary);
  }

/*
  // Get some specific element in the mesh and check if it is on the boundary of not
  printf("checking if some entity is on the boundary\n");
  int face_number = 0;
  pMeshEnt face_tocheck = pumi_mesh_findEnt(mesh,2,face_number);
  // get the adjacent edges of the face
  Adjacent adjacent;
  int num_adj = pumi_ment_getAdjacent(face_tocheck,1,adjacent);
  for (int i = 0;i<num_adj;i++){
    int boundary_number = get_edge_bound_num(adjacent[i],all_on_boundary);
    if(boundary_number){
      printf("Is on boundary %d\n", boundary_number);
      break;
    }
    else if(i == num_adj-1){
      printf("Is not on boundary \n");
    }
  }

  // Check if the get face area routine is working
  double myarea = get_face_area(face_tocheck);
  printf("The area is %f \n", myarea);
  pMeshEnt e;
  double totalarea = 0;
  pMeshIter it = mesh->begin(2);
  while ((e = mesh->iterate(it))){
    totalarea = totalarea + get_face_area(e);
  }
  mesh->end(it);
  printf("Total area is %f \n", totalarea);
*/



  // Use some library to solve the assembled matrix


  // Now run some code to generate the secondary variables and assign them to some field


  // Visualize the solution with Paraview directly if possible




  pumi_finalize();
  MPI_Finalize();
  cout << "End of Program\n";
}
