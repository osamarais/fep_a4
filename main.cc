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
  get_element_type(mesh);

  // Use some library to solve the assembled matrix


  // Now run some code to generate the secondary variables and assign them to some field


  // Visualize the solution with Paraview directly if possible




  pumi_finalize();
  MPI_Finalize();
  cout << "End of Program\n";
}
