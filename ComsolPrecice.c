#include "precice/SolverInterfaceC.h"
#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>

#include "fsi_mesh.h"
#include "fsi_mesh_utils.h"
#include "fsi_map-mat.h"
#include "fsi_interface_socket.h"

// on this socket is the communication between COMSOL and preCICE trough FSI*ce
#define COMMUNICATION_SOCKET       53219
#define HOSTNAME                  "localhost"
#define MAX_NR_POINTS              10000

static int socketToComsol;
static FSI_Mesh * mesh = NULL;
static FSI_Data * meshForces = NULL;
static FSI_Data * meshDisplacements = NULL;
static FSI_Data * meshDisplacementDeltas = NULL;
static FSI_Data * meshVelocities = NULL;
static FSI_Data * meshVelocityDeltas = NULL;

static int iterationNumber;
static int maxIterationNumber = 20;
static double timeStep = 0.01;
static int redoStep = 0;
static int doStep = 1;
static int succesStep = 1;
static int forcesID = -1;
static int velocitiesID = -1;
static int displacementsID = -1;
static int displacementDeltasID = -1;
static int velocityDeltasID = -1;
static int meshID = -1;
static int dimensions = 0;

//static double x_actualshift[MAX_NR_POINTS];
//static double y_actualshift[MAX_NR_POINTS];
static FILE *fd;


void readForces
(
  FSI_Mesh * mesh,
  FSI_Data * data,
  int        dataID )
{
  int i;
  double coords[dimensions];
  double value[dimensions];

  fprintf ( fd, "Reading forces (dim = %d)\n", dimensions );

  if ( dimensions == 2 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh, i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      //@todo replace hard-coded vertexIDs "i-1" with the ones received from the mesh definition
      precicec_readVectorData ( dataID, i-1, value );
      FSI_Data_set_vector ( data, i, value[0], value[1], 0.0 );
    }
  }
  else if ( dimensions == 3 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh , i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      coords[2] = FSI_Mesh_get_node_z ( mesh, i );
      //@todo replace hard-coded vertexIDs "i-1" with the ones received from the mesh definition
      precicec_readVectorData ( dataID, i-1, value );
      FSI_Data_set_vector ( data, i, value[0], value[1], value[2] );
    }
  }
}

void writePreciceData
(
  FSI_Mesh * mesh,
  FSI_Data * data,
  int        dataID )
{
  int i;
  int iDim;
  double coords[dimensions];
  double value[dimensions];

  fprintf ( fd, "Writing data \"%s\"\n", data->_name );

  if ( dimensions == 2 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh, i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      value[0] = FSI_Data_get_value ( data, i, 1 );
      value[1] = FSI_Data_get_value ( data, i, 2 );
      fprintf ( fd, "Set value %d = %e,%e at coords = %e, %e\n",
                i, value[0], value[1], coords[0], coords[1] );
      //@todo replace hard-coded vertexIDs "i-1" with the ones received from the mesh definition          
      precicec_writeVectorData ( dataID, i-1, value );
    }
  }
  else if ( dimensions == 3 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh, i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      coords[2] = FSI_Mesh_get_node_z ( mesh, i );
      value[0] = FSI_Data_get_value ( data, i, 1 );
      value[1] = FSI_Data_get_value ( data, i, 2 );
      value[2] = FSI_Data_get_value ( data, i, 3 );
      fprintf ( fd, "Set value %d = %e,%e, %e at coords = %e, %e, %e\n",
                i, value[0], value[1], value[2], coords[0], coords[1], coords[2] );
      //@todo replace hard-coded vertexIDs "i-1" with the ones received from the mesh definition
      precicec_writeVectorData ( dataID, i-1, value );
    }
  }
}


//====================================== MAIN =========================================
int main(int argc, char** argv)
{
  int hServerSocket, hServerSocket1; /* handle to socket */
  struct hostent * pHostInfo;        /* holds info about a machine */
  struct sockaddr_in Address;       /* Internet socket address stuct */
  long nHostAddress;       /* long address of the host we want to connect to */
  int fHostPort, sHostPort;
  int nAddressSize = sizeof(struct sockaddr_in);
	double coords[2], value[2];
//	int pointNr , triangleNr;
	int i;
  int checkpointNumber;
	char filename[50];

	fd = fopen("comsolprecice.tmp", "w");

  if (argc < 1){
    // fprintf(fd,"TOO FEW PARAMETERS: \n -accessorName [IN] Name of the solver accessing the interface. Has to match one of the names specified in the configuration xml file\n");
    fprintf(fd,"TOO FEW PARAMETERS: \n");
    fprintf(fd," -configFileName [IN] (Path and) name of the xml configuration file containing the precice configuration.");
	  return 0;
	}

  fprintf(fd, "Calling precicec_createSolverInterface with: Comsol, %s \n", argv[1]);
  precicec_createSolverInterface("Comsol", argv[1], 0, 1);
  dimensions = precicec_getDimensions();

  socketToComsol = socket(AF_INET, SOCK_STREAM, 0);
  pHostInfo = gethostbyname(HOSTNAME); /* get IP address from name */
  /* copy address into long */
  memcpy(&nHostAddress, pHostInfo->h_addr, pHostInfo->h_length);
  /* fill address struct */
  Address.sin_addr.s_addr = nHostAddress;
  Address.sin_port = htons(COMMUNICATION_SOCKET);
  Address.sin_family = AF_INET; /* AF_INET represents the address family INET for Internet sockets. */
  fprintf(fd,"Connecting to %d on port %d \n", nHostAddress, sHostPort);
  /* connect to host */
  printf("Connecting to Comsol...\n"); fflush(stdout);
  while(connect(socketToComsol, (struct sockaddr*) &Address, sizeof(Address)) == -1){
    usleep(1000000);
  }
	printf("Connection with success ... Now receiving Mesh\n"); fflush(stdout);

  // --- send the checkpoint number which should be loaded ---
  checkpointNumber = -1; // -1 means it no checkpoints will be loaded
                          // only when checkpoint_number is greater than 0,
                          // will a checkpoint loaded
  FSI_Send_int_through_socket(socketToComsol, &checkpointNumber);

  /* GET THE MESH FROM COMSOL AND MAKE INITIALIZATIONS */
	mesh = FSI_Mesh_new_empty(); /* this initialization is needed */
	fprintf(fd, "Receiving FSIce mesh ... \n");
	FSI_Recv_mesh_socket(mesh , socketToComsol);

  fprintf(fd, "Building preCICE mesh ... \n");
  int meshID = precicec_getMeshID("Comsol-Mesh");

	/* create the vertices */
	for (i=1; i < FSI_Mesh_get_num_nodes(mesh); i++){
    coords[0] = FSI_Mesh_get_node_x(mesh, i);
    coords[1] = FSI_Mesh_get_node_y(mesh, i);
    fprintf(fd, "Adding vertex with coodrs = %f, %f\n", coords[0], coords[1]);
    //@todo store index
    int index = precicec_setMeshVertex(meshID, coords);
    fprintf(fd, "   ... returned index = %d\n", index);
	}
	/* Create the edges */
	for (i=0; i < FSI_Mesh_get_num_triangles(mesh); i++){
	  fprintf(fd, "Adding edge with vertex indices = %d, %d\n",
	          FSI_Mesh_get_face_node1(mesh,i)-1, FSI_Mesh_get_face_node2(mesh,i)-1);
    //@todo replace hard-coded vertexIDs with stored ones
	  precicec_setMeshEdge(meshID,
                        FSI_Mesh_get_face_node1(mesh,i)-1,
                        FSI_Mesh_get_face_node2(mesh,i)-1);
	}

  fprintf(fd, "Receive and create data ...\n");
  meshForces = FSI_Mesh_Data_new(mesh, 1, 3, "Forces");
	FSI_Recv_quantity_socket(mesh, socketToComsol);
  FSI_Recv_quantity_socket(mesh, socketToComsol);
  FSI_Recv_quantity_socket(mesh, socketToComsol);
  FSI_Recv_quantity_socket(mesh, socketToComsol);
	meshDisplacements = FSI_Mesh_get_Data(mesh, "Displacements");
  meshDisplacementDeltas = FSI_Mesh_get_Data(mesh, "DisplacementDeltas");
  meshVelocities = FSI_Mesh_get_Data(mesh, "Velocities");
  meshVelocityDeltas = FSI_Mesh_get_Data(mesh, "VelocityDeltas");

  printf("Initializing Coupling ... \n"); fflush(stdout);
	timeStep = precicec_initialize();

	forcesID = precicec_getDataID("Forces");
	if (precicec_hasData("Velocities", meshID)){
	  velocitiesID = precicec_getDataID("Velocities", meshID);
	}
	if (precicec_hasData("VelocityDeltas", meshID)){
	  velocityDeltasID = precicec_getDataID("VelocityDeltas", meshID);
	}
	if (precicec_hasData("Displacements", meshID)){
	  displacementsID = precicec_getDataID("Displacements", meshID);
	}
	if (precicec_hasData("DisplacementDeltas", meshID)){
	  displacementDeltasID = precicec_getDataID("DisplacementDeltas", meshID);
	}

  // MAIN TIMESTEPPING LOOP
  //
  // Overview:
  //
  // if (do)
  //    send do to comsol
  //    send redo to comsol
  //    send dt to comsol
  //    copy forces from precice to fsice mesh
  //    send forces to comsol
  //    receive comsol values
  //    copy comsol values from fsice to precice mesh
  //    send do checkpoint to comsol
  // send do to comsol
  // send do checkpoint to comsol
  //
  while (precicec_isCouplingOngoing() > 0){
    iterationNumber++;
	  if (precicec_isActionRequired("write-iteration-checkpoint")){
      precicec_markActionFulfilled("write-iteration-checkpoint");
	  }
	  if (precicec_isActionRequired("read-iteration-checkpoint")){
		  redoStep = 1; /* repeat the timestep */
	    precicec_markActionFulfilled("read-iteration-checkpoint");
	  }
	  else {
	    redoStep = 0; /* do not repeat the timestep */
	  }
	  fprintf(fd, "Sending coupling state to Comsol ...\n");
    FSI_Send_int_through_socket(socketToComsol, &doStep);
    FSI_Send_int_through_socket(socketToComsol, &redoStep);
    FSI_Send_double_through_socket(socketToComsol, &timeStep);

    readForces(mesh, meshForces, forcesID);
    fprintf(fd, "Sending force data to Comsol ...\n");
    FSI_Send_quantity_socket(mesh, "Forces", socketToComsol);

	  fprintf(fd, "Receiving Comsol data 1 ...\n");
	  FSI_Recv_int_through_socket(socketToComsol , &succesStep);
	  fprintf(fd, "Receiving Comsol data 2 ...\n");
	  FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "Receiving Comsol data 3 ...\n");
    FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "Receiving Comsol data 4 ...\n");
    FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "Receiving Comsol data 5 ...\n");
    FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "... done\n");

    if (velocitiesID != -1){
      writePreciceData(mesh, meshVelocities, velocitiesID);
    }
    if (velocityDeltasID != -1){
      writePreciceData(mesh, meshVelocityDeltas, velocityDeltasID);
    }
    if (displacementsID != -1){
	    writePreciceData(mesh, meshDisplacements, displacementsID);
    }
    if (displacementDeltasID != -1){
	    writePreciceData(mesh, meshDisplacementDeltas, displacementDeltasID);
    }

	  // ------- Exchange data -----------
    fprintf ( fd, "Calling preCICE advance ... \n");
	  timeStep = precicec_advance(timeStep);
	  //@todo extend to sub-cycling

  }
	fprintf(fd,"Tell COMSOL simulation ended  ... \n");

	doStep = 0;
  FSI_Send_int_through_socket ( socketToComsol , &doStep );
  close ( socketToComsol );
	fprintf ( fd, "End simulation ... \n" );
	fclose ( fd );
  return 1;
}
