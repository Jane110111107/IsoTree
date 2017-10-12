#include "LingoProgram.h"
#include <stdio.h>
#include <string.h>

/*
int CALLTYPE  MyCallback( void* pModel, int nReserved, void* pUserData) {

   // this callback function will be called periodically by the Lingo solver 

   int* pnCallbacks = (int*) pUserData;
   ++*pnCallbacks;
   return( 0);

}
*/


void lingo_solve1() {
printf("lingo1 begin\n");
   char cScript[256];
   int nError=-1;//, nPointersNow;
  // int nCallbacks = 0;
  // double  dStatus=-1.;

   // create the LINGO environment object
   pLSenvLINGO pLINGO;
   pLINGO = LScreateEnvLng();

   if ( !pLINGO) {

      printf( "Can''t create LINGO environment!\n");
      goto FinalExit;

   }

   // Pass LINGO a pointer to our callback function
  // nCallbacks = 0;
  // nError = LSsetCallbackSolverLng( pLINGO, &MyCallback, &nCallbacks);
  // if ( nError) goto ErrorExit;

   // Pass LINGO a pointer to our callback function
   //nCallbacks = 0;
   //nError = LSsetCallbackErrorLng( pLINGO, &MyErrorCallback, NULL);
  // if ( nError) goto ErrorExit;

     // Open LINGO's log file  
   //nError = LSopenLogFileLng( pLINGO, "LINGO.log");
   //if ( nError) goto ErrorExit;

 // nError = LSsetPointerLng( pLINGO, &dStatus, &nPointersNow); 
 // if ( nError) goto ErrorExit;
 // printf("strcpy\n");
   // Here is the script we want LINGO to run
   strcpy( cScript, "SET ECHOIN 1 \n TAKE mip1.lng \n GO \n QUIT \n");
  printf("strcpy end\n");
   // Run the script
   nError = LSexecuteScriptLng( pLINGO, cScript);
   if ( nError) goto ErrorExit;

   // Close the log file
   //LScloseLogFileLng( pLINGO);

   // Any problems?
   if ( nError )
	   printf( "\nUnable to solve!!!\n\n");

   goto NormalExit;

ErrorExit:
   printf("LINGO Error Code: %d\n", nError);

NormalExit:
   LSdeleteEnvLng( pLINGO);

FinalExit: ;
printf("lingo1 end\n");
}
void lingo_solve2() {
printf("lingo1 begin\n");
   char cScript[256];
   int nError=-1;//, nPointersNow;
  // int nCallbacks = 0;
   //double  dStatus=-1.;

   // create the LINGO environment object
   pLSenvLINGO pLINGO;
   pLINGO = LScreateEnvLng();

   if ( !pLINGO) {

      printf( "Can''t create LINGO environment!\n");
      goto FinalExit;

   }

   // Pass LINGO a pointer to our callback function
   //nCallbacks = 0;
 //  nError = LSsetCallbackSolverLng( pLINGO, &MyCallback, &nCallbacks);
  // if ( nError) goto ErrorExit;

   // Pass LINGO a pointer to our callback function
   //nCallbacks = 0;
   //nError = LSsetCallbackErrorLng( pLINGO, &MyErrorCallback, NULL);
  // if ( nError) goto ErrorExit;

     // Open LINGO's log file  
   //nError = LSopenLogFileLng( pLINGO, "LINGO.log");
   //if ( nError) goto ErrorExit;

  //nError = LSsetPointerLng( pLINGO, &dStatus, &nPointersNow); 
  //if ( nError) goto ErrorExit;
  
   // Here is the script we want LINGO to run
   strcpy( cScript, "SET ECHOIN 1 \n TAKE mip2.lng \n GO \n QUIT \n");

   // Run the script
   nError = LSexecuteScriptLng( pLINGO, cScript);
   if ( nError) goto ErrorExit;

   // Close the log file
   //LScloseLogFileLng( pLINGO);

   // Any problems?
   if ( nError )
	   printf( "\nUnable to solve!!!\n\n");

   goto NormalExit;

ErrorExit:
   printf("LINGO Error Code: %d\n", nError);

NormalExit:
   LSdeleteEnvLng( pLINGO);

FinalExit: ;
printf("lingo2 end\n");
}



