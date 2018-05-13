/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat    
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>



int main(void)
{  

    // Initialisation des grains
    int n = 45;
    double radius    = 0.05;
    double mass      = 0.1;
    double radiusIn  = 0.5;
    double radiusOut = 2.0;    
    double dt      = 1e-2;
    double tEnd    = 8.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100;
    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut);

    // Initialisation de la vitesse d'exécution des mouvements des grains
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.005;

    // Solveur bande en numérotation selon x
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_XNUM;

    // Chargement du maillage
    char meshFileName[] = "../data/meca1120-projet-meshMedium.txt";   
    
    // Initialisation du problème
    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);


    clock_t tic = clock();

    // Résolution du problème
    femDiffusionCompute(theProblem);

    // Print des infos sur le maillage - inutile dans le projet
    femSolverPrintInfos(theProblem->solver);


    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    

    int option = 1;
    // Récupération si changement de type de solveur
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("MECA1120 : Projet 2018");
    glfwMakeContextCurrent(window);

    do 
    {
        int testConvergence,w,h,i;
        double currentTime = glfwGetTime();
        char theMessage1[256],theMessage2[256];
        sprintf(theMessage1, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        sprintf(theMessage2, "Time : %g sec",t);
        glfwGetFramebufferSize(window,&w,&h);
      
        if (option == 1) {
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotField(theProblem->mesh,theProblem->soluce);   
        }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); 
        }

        for (i=0 ;i < theGrains->n; i++) {     
            glColor3f(1,0,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); 
        }       
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);

        glColor3f(0.0,0.0,0.0); glfemDrawMessage(20,460,theMessage1);
        glColor3f(0.0,0.0,0.0); glfemDrawMessage(100,460,theMessage2);

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);  
  //
  // A decommenter pour pouvoir progresser pas par pas
  //          printf("press CR to compute the next time step >>");
  //          char c= getchar();
  //
            femGrainsUpdate(theGrains,dt,tol,iterMax);
            t += dt; 
        }

        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS) 
                theRunningMode = 1; 
          if (glfwGetKey(window,'S') == GLFW_PRESS) 
                theRunningMode = 0; 
        }   
      
      /*
        if (solverType != newSolverType || renumType != newRenumType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
    */

    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );

            
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
}

