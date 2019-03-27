//#include "reorder.h"
#include <queue>
// This file contains functions to reorder the mesh, with headers in the header file





pMeshEnt getStart(pMesh &mesh, pGeom &g){
  pGeomEnt mystartg;
  int minedgesyet = 0;
  pMeshEnt mystarte;
  for (pGeomIter it = g->begin(0); it!=g->end(0);++it){
    pGeomEnt e = *it;

    // Fill in on first iteration
    if (0 == minedgesyet){
      mystartg = e;
      minedgesyet = pumi_gent_getNumAdj(e,1);
    }
    // Check on next iterations
    if(minedgesyet > pumi_gent_getNumAdj(e,1) ){
      mystartg = e;
      minedgesyet = pumi_gent_getNumAdj(e,1);
    }
    // For a 2D case, 2 vertices is the minimum. Exit if this happens.
    if(2 == minedgesyet){
      break;
    }
  }
  // Get reverse classification of mesh entity on geometic vertex
  // return the mesh vertex
  std::vector<pMeshEnt> myreversevector;
  pumi_gent_getRevClas(mystartg,myreversevector);
  return myreversevector[0];
  //printf("Got the starting vertex\n");
}












void reorder(pMesh &mesh, pMeshEnt &startingvertex){
  pShape myshape = pumi_mesh_getShape(mesh);
  pNumbering mynum =  pumi_numbering_create(mesh,"mynewnumbering",myshape);
  //pMeshTag facenum =  pumi_mesh_createIntTag(mesh,"myfacenumbering",1);
  pShape myfaceshape = pumi_shape_getConstant(2);
  pNumbering facenum = pumi_numbering_create(mesh,"myfacenumbering",myfaceshape);

  // First get the total number of nodes
  int labelnode = 1;
  labelnode += (pumi_shape_getNumNode(myshape, 0) * pumi_mesh_getNumEnt(mesh,0));
  labelnode += (pumi_shape_getNumNode(myshape, 1) * pumi_mesh_getNumEnt(mesh,1));
  //printf("The total number of nodes is %d\n", labelnode-1);
  // Now get the total number of faces since we will be labelling them also
  int labelface = pumi_mesh_getNumEnt(mesh,2) + 1;
  //printf("The total number of faces is %d\n", labelface-1);
  // Create the queue and the list
  std::deque<pMeshEnt> myq;
  std::deque<pMeshEnt> mylist;
  // add the starting vertex to the queue
  myq.push_back(startingvertex);



  while(myq.size() > 0){
    pMeshEnt entity = myq.front();
    //myq.pop_front();
    // If node associated with this mesh entity is unlabelled then label it
    if (!pumi_node_isNumbered(mynum,entity,0,0)){
      labelnode = labelnode - 1;
      pumi_node_setNumber(mynum,entity,0,0,labelnode);
      //printf("LABELLED A NODE  %d\n", labelnode);
    }
    // If this entity is a vertex
    // Find the edges adjacent to the vertex
    // Find the faces adjacent to the edges
    // number these faces
    if(0==pumi_ment_getDim(entity)){
      // get the adjacent edges
      Adjacent myadje;
      int numedges = pumi_ment_getAdjacent(entity,1,myadje);
      //printf("The number of adjacent edges is %d \n", numedges);
      // loop over these adges
      for (int i = 0; i < numedges; i++){
        //printf("looping over adjacent edges %d\n", i);
        pMeshEnt thisedge = myadje[i];
        // get the adjacent faces
        Adjacent myadjf;
        int numfaces = pumi_ment_getAdjacent(thisedge,2,myadjf);
        //printf("The number of adjaent faces is %d \n", numfaces);
        for (int j = 0; j < numfaces; j++){
          //printf("looping over faces %d \n", j);
          if (!pumi_node_isNumbered(facenum,myadjf[j],0,0)){
            labelface = labelface - 1;
            pumi_node_setNumber(facenum,myadjf[j],0,0,labelface);
            //printf("LABELLED FACE %d\n", labelface);
          }
          // If face has node that needs queueing queue it
          // There is no such case in our given meshes, simply skip this part
          // Implement this later on
        }



        // othervertex
        pMeshEnt othervertex = pumi_medge_getOtherVtx(thisedge,entity);
        //printf("Got other vertex\n");
        if (pumi_shape_hasNode(myshape,1)){
          // if othervertex labelled or in queue and edge node not labelled
          // check if othervertex is in queue
          bool isinq;
          std::deque<pMeshEnt>::iterator it = myq.begin();
          while(it != myq.end()){
            if(*it == othervertex){
              isinq = true;
              break;
            }
            else {
              isinq = false;
            }
            it++;
          }
          bool isotherlabelled = pumi_node_isNumbered(mynum,othervertex);
          bool isedgelabelled = pumi_node_isNumbered(mynum,thisedge);
          if ((isinq or isotherlabelled) and !isedgelabelled){
            labelnode = labelnode - 1;
            pumi_node_setNumber(mynum,thisedge,0,0,labelnode);
            //printf("LABELLED NODE %d\n", labelnode);
          }
          else{
            if (!isedgelabelled){
              myq.push_back(thisedge);
            }
            if(!isotherlabelled){
              mylist.push_back(othervertex);
            }
          }
        }
        else{
          bool isotherlabelled = pumi_node_isNumbered(mynum,othervertex);
          if(!isotherlabelled){
            //printf("The other vertex was not labeller, so put in mylist\n");
            mylist.push_back(othervertex);
          }
        }


        // Closed Loop Over Edges
      }
    }

    // now empty the list into the queue
    //printf("The size of the list is %d \n", mylist.size());
    while (mylist.size() != 0){
      //printf("Looping in list\n");
      pMeshEnt check = mylist.front();
      mylist.pop_front();
      bool isinq;
      // add elements into the queue that are already not in there
      std::deque<pMeshEnt>::iterator it = myq.begin();
      while(it != myq.end()){
        //printf("Iteratinf in the queue\n");
        if(*it == check){
          isinq = true;
          //printf("Found it in the list \n");
          break;
        }
        else{
          //printf("Did not find it in the list\n");
          isinq = false;
        }
        it++;
      }
      if(isinq == false){
        myq.push_back(check);
        //printf("Added an element to the queue\n");
      }
    }
    myq.pop_front();
  }
}



void reorder_mesh(pMesh mesh, pGeom geom){
  pMeshEnt startingvertex = getStart(mesh,geom);
  reorder(mesh,startingvertex);
}





bool hasNode(pMesh m, pMeshEnt e)
{
  return pumi_shape_hasNode(pumi_mesh_getShape(m),pumi_ment_getTopo(e));
}
