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
  // Try to find which boundary something is on
  std::vector<int> list;
  int tofind = 0;
  pMeshEnt ment = pumi_mesh_findEnt(mesh,0,tofind);
  get_bound_num(ment, boundary_verts, list);
  for (int i = 0; i< list.size(); i++){
    printf("The vertex %d is on %d\n",tofind,list[i]);
  }





  // Use some library to solve the assembled matrix


  // Now run some code to generate the secondary variables and assign them to some field


  // Visualize the solution with Paraview directly if possible




  pumi_finalize();
  MPI_Finalize();
  cout << "End of Program\n";
}
