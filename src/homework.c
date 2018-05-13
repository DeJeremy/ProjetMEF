// Devoir 5 réalisé par Degraeuwe Jérémy
// NOMA : 23891600
// Attention : l'idée de regrouper les coordonnées et l'indice d'un noeud au sein d'une structure a été trouvée en collaboration avec Jimmy Fraiture

#include"fem.h"

typedef struct 
{
    double X;
    double Y;
    int indice;
} renum;

void renumCreate(renum *r,femMesh *mesh,int n)
{
    int i;
    for(i=0 ; i<n ; i++)
    {
        r[i].X = mesh->X[i];
        r[i].Y = mesh->Y[i];
        r[i].indice = i;
    }
}

int elemCompareX(const void *e1, const void *e2)
{
    if(((renum*)e1)->X < ((renum*)e2)->X)
        return 1;
    else if(((renum*)e1)->X > ((renum*)e2)->X)
        return -1;
    else
        return 0;
}

int elemCompareY(const void *e1, const void *e2)
{
    if(((renum*)e1)->Y < ((renum*)e2)->Y)
        return 1;
    else if(((renum*)e1)->Y > ((renum*)e2)->Y)
        return -1;
    else
        return 0;
}

int ind(int *number, int oldInd, int n)
{
    int i;
    for(i=0 ; i<n ; i++)
    {
        if(number[i] == oldInd)
            return i;
    }
    return 0;
}


void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i;
    renum *elemsRenum = malloc(sizeof(renum)*theProblem->mesh->nNode);
    renumCreate(elemsRenum,theProblem->mesh,theProblem->mesh->nNode);
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
            break;
        case FEM_XNUM : 
            qsort(elemsRenum, theProblem->mesh->nNode, sizeof(renum),elemCompareX);
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[elemsRenum[i].indice] = i;
            break;
        case FEM_YNUM : 
            qsort(elemsRenum, theProblem->mesh->nNode, sizeof(renum),elemCompareY);
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[elemsRenum[i].indice] = i;
            break;          

        default : Error("Unexpected renumbering option"); }

        free(elemsRenum);
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    int myBand      = 0,
        nElem       = theMesh->nElem,
        nLoc        = theMesh->nLocalNode,
        currentBand = 0,
        i,
        j,
        k,
        l,
        m;
        for(k=0 ; k<nElem ; k++)
        {
            for(m=0 ; m<nLoc ; m++)
            {
                for(l=m+1 ; l<nLoc ; l++)
                {
                    i = theProblem->number[theMesh->elem[(k*nLoc)+m]];
                    j = theProblem->number[theMesh->elem[(k*nLoc)+l]];
                    currentBand = abs(i-j);
                    if(currentBand > myBand)
                        myBand = currentBand;
                }
            }
        }
    return myBand+1;
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    double *R   = mySolver->R,
           *D   = mySolver->D,
           *S   = mySolver->S;
    int i,
        j,
        myRow,
        myRow2;

    if(mySolver->iter == 0)
    {
        for (i=0 ; i<nLoc ; i++) 
        { 
            myRow = map[i];
            for(j=0 ; j<nLoc ; j++) 
            {
                R[myRow] += Aloc[(i * nLoc) + j] * Uloc[j];
            }
            R[myRow] += - Bloc[i];
        }
    }

    for(i=0 ; i<nLoc ; i++)
    {
        myRow = map[i];
        for(j=0 ; j<nLoc ; j++)
        {
            myRow2 = map[j];
            S[myRow] += Aloc[(i * nLoc) + j] * D[myRow2];
        }
    }
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) 
{   
    mySolver->R[myNode] = 0.0;
    mySolver->S[myNode] = 0.0;
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    double error    = 0.0,
           alpha    = 0.0,
           beta     = 0.0,
           *R       = mySolver->R,
           *D       = mySolver->D,
           *S       = mySolver->S,
           *X       = mySolver->X,
           size     = mySolver->size;
    int i;

    for (i=0 ; i<size ; i++)
        error += R[i] * R[i];

    if(mySolver->iter == 0)
    {
        for (i=0 ; i<size ; i++)
        {
            X[i] = 0.0;
            D[i] = R[i];
        }
    }
    else
    {
        for (i=0 ; i<size ; i++)
            alpha += R[i] * S[i];
        alpha = - error / alpha;

        for (i=0 ; i<size ; i++)
        {
            R[i] += alpha * S[i];
            beta += R[i] * R[i];
        }
        beta = beta / error;

        for (i=0 ; i<size ; i++)
        {
            X[i] = alpha * D[i];
            D[i] = R[i] + (beta * D[i]);
            S[i] = 0.0;
        }
    }

    mySolver->iter++;
    mySolver->error = sqrt(error);
    return(X);
}

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  
{
    int     n           = myGrains->n,
            i,
            j,
            iContact;
    double  *x          = myGrains->x,
            *y          = myGrains->y,
            *m          = myGrains->m,
            *r          = myGrains->r,
            *vy         = myGrains->vy,
            *vx         = myGrains->vx,
            *dvBoundary = myGrains->dvBoundary,
            *dvContacts = myGrains->dvContacts,
            rIn         = myGrains->radiusIn,
            rOut        = myGrains->radiusOut,
            error       = 0.0,
            gap,rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;

    error = 0;
    iContact = 0;
    for(i=0 ; i<n ; i++) 
    {
        for(j=i+1 ; j<n ; j++) 
        { 
            rx = (x[j]-x[i]);
            ry = (y[j]-y[i]);
            rr = sqrt(rx*rx+ry*ry);
            nx = rx/rr;
            ny = ry/rr;
            if (iter == 0) 
            {
                  dv = dvContacts[iContact]; 
            }
            else 
            {
                  vn = (vx[i]-vx[j])*nx + (vy[i]-vy[j])*ny ;
                  gap = rr - (r[i]+r[j]);
                  dv = fmax(0.0, vn + dvContacts[iContact] - gap/dt);
                  dv = dv - dvContacts[iContact];                      
                  dvContacts[iContact] += dv; 
                  error = fmax(fabs(dv),error); 
            }
            vx[i] -= dv * nx * m[j] / ( m[i] + m[j] );
            vy[i] -= dv * ny * m[j] / ( m[i] + m[j] );
            vx[j] += dv * nx * m[i] / ( m[i] + m[j] );
            vy[j] += dv * ny * m[i] / ( m[i] + m[j] );                  
            iContact++; 
        }
    }

    for(i=0 ; i<n ; i++) 
    {
        rr = sqrt(x[i] * x[i] + y[i] * y[i]);      
        nx = x[i] / rr;
        ny = y[i] / rr;
        if (iter == 0) 
        {
            dv = dvBoundary[i]; 
        }
        else 
        {
            vn = vx[i]*nx + vy[i]*ny ;      
            gap = rOut - rr - r[i];
            dvOut = fmax(0.0, vn + dvBoundary[i] - gap/dt);
            gap = rr - rIn - r[i];
            dvIn  = fmax(0.0,- vn - dvBoundary[i] - gap/dt);
            dv = dvOut - dvIn - dvBoundary[i]; 
            dvBoundary[i] += dv; 
            error = fmax(fabs(dv),error);
        }
        vx[i] -= dv * nx;
        vy[i] -= dv * ny;
    }  
    return error;
}

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax)
{
    int     i,
            n           = myGrains->n,
            iter        = 0;
    double  *x          = myGrains->x,
            *y          = myGrains->y,
            *m          = myGrains->m,
            *vy         = myGrains->vy,
            *vx         = myGrains->vx,
            gamma       = myGrains->gamma,
            gx          = myGrains->gravity[0],
            gy          = myGrains->gravity[1],
            error;

// 
// -1- Calcul des nouvelles vitesses des grains sur base de la gravite et de la trainee
//

    for(i=0 ; i<n ; i++) 
    {
        //double fx = m[i] * gx - gamma * vx[i]; // fx = m[i] * gx - (gamma * (vx[i] - u(x[i],y[i]))
        //double fy = m[i] * gy - gamma * vy[i]; // fy = m[i] * gy - (gamma * (vy[i] - v(x[i],y[i]))
        // ========== tests ==========
        double fx = m[i] * gx - gamma * vx[i];
        double fy = m[i] * gy - gamma * vy[i];
        // ==========      ===========
        vx[i] += fx * dt / m[i];
        vy[i] += fy * dt / m[i];  
    }

//
// -2- Correction des vitesses pour tenir compte des contacts        
//       
           
    do 
    {
        error = femGrainsContactIterate(myGrains,dt,iter);
        iter++; 
    }
    while ((error > tol/dt && iter < iterMax ) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,error);
 
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//

    for (i=0 ; i<n ; i++) 
    {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; 
    }
}

