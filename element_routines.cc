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
}



// The ordering of the elements must be ascertained properly to assemble the matrix.
// This means that: loop over the elements in any order that is found to be okay
// Simply get the numbering from the appropriate numbering we have created for the faces and the nodes.
// Over here we must check if we have already done the reordering or not so that we may use the appropriate numbering from the mesh


// Write element routines here to get the coefficients of the elemental matrices.



// Write some functions here to deal with the boundary conditions.



// Write a function to calculate the areas of the element if required
