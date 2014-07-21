/*
 * Copyright (c) David Powell <david@drp.id.au>
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 */


// Ukkonen, linear, checkpointing
//
//
//   Date: 12/9/96
//   Author: David Powell

#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>

#include "common.h"
#include "ukk_linear.h"

using namespace std;

#undef MODSIZE
#define MODSIZE      (MAX2((a+b)+1, MismatchCost)+1)


#define MAXDIR3(x,y,z) ((x)>=(y) ? \
			((x)>=(z) ? diag : vert) : \
			((y)>=(z) ? horz : vert))


#define CHECKPSIZE   (MAX2((a+b), MismatchCost)+1)
#define NUMDIRS 3		// Three directions horz,vert,diag.

class UkkonenCheckp
{
  struct ukkElem {
    int val;		     // Contents of this cell
    int cost;		     // Cost of this cell (used for debugging)
    int done;		     // Flag for whether this cell has been calculated
  };

  enum direction {nodir=-1, diag, horz, vert};

  struct checkpElem {
    int cost, diag;
    int dist;
    
    direction entryDir;		// Direction moved into this cell

    direction fromDir;
    int fromCost, fromDiag;    
  };
  
  
  char *A, *B;			// Two strings to be aligned
  int lenA,lenB;		// Length of the two strings
  struct ukkElem (*data)[MODSIZE][NUMDIRS];

  struct checkpElem (*checkp)[CHECKPSIZE][NUMDIRS];
  int checkpCost;
  
  int offset;			// Offset into data[] because -ve indices
				// needed.
  int arraySize;		// Number of rows allocated for each array

#if DEBUG
  int depth=-2;
#endif

private:
  void storeCheckp(int c, int d, direction dir, int dist,
		   int fromCost, int fromDiag, direction fromDir)
  {
    checkp[d+offset][c%CHECKPSIZE][dir].cost = c;
    checkp[d+offset][c%CHECKPSIZE][dir].diag = d;
    checkp[d+offset][c%CHECKPSIZE][dir].entryDir = dir;
    checkp[d+offset][c%CHECKPSIZE][dir].dist = dist;
    checkp[d+offset][c%CHECKPSIZE][dir].fromCost = fromCost;
    checkp[d+offset][c%CHECKPSIZE][dir].fromDiag = fromDiag;
    checkp[d+offset][c%CHECKPSIZE][dir].fromDir  = fromDir;
  }

  int Ukk(int sDiag, int sCost, direction matrix, int d, int c)
  {
    int v1,v2,v3,res;
    int fromCost, fromDiag;
    direction fromDir;
  
//    if (d==sDiag && c==sCost) return sDistance; // Starting condition

#if DEBUG
    for (v1=0;v1<depth+2;v1++)
      printf(" ");
    printf("D=%d C=%d DIR=%s\n",d,c,
	   (matrix==diag ? "Diag" : (matrix==horz ? "Horz" : "Vert")));
#endif
    
    
    if (abs(d-sDiag) > c-sCost)
      return BIG_NEGATIVE;	//Boundary condition

    if (data[d+offset][c%MODSIZE][matrix].done)	// Check if already calculated
      if (data[d+offset][c%MODSIZE][matrix].cost==c)
        return data[d+offset][c%MODSIZE][matrix].val;
      else {
        fprintf(stderr,"INTERNAL ERROR: Use of column that is not correct");
        fprintf(stderr,"\nCost = %d,  Diag = %d, (sCost=%d)  (Cost in Matrix=%d)\n",
                c,d,sCost,data[d+offset][c%MODSIZE][matrix].cost);
        exit(-1);
      }

#ifdef DEBUG    
    depth+=2;
    for (v1=0;v1<depth;v1++)
      printf(" ");
    printf("Depends on:\n");
#endif
    
    // Must calculate if not previously done.
    switch (matrix) {
     case diag : 
      v1 = Ukk(sDiag,sCost, diag, d, c-MismatchCost)+1; 
      v2 = Ukk(sDiag,sCost, horz, d, c-MismatchCost)+1; 
      v3 = Ukk(sDiag,sCost, vert, d, c-MismatchCost)+1;
      switch (fromDir=MAXDIR3(v1,v2,v3)) {
      case diag : res = v1; fromCost = c-MismatchCost; break;
      case horz : res = v2; fromCost = c-MismatchCost; break;
      case vert : res = v3; fromCost = c-MismatchCost; break;
      default :	assert(0);
      }

      // Need to also check same cost for horz followed by a run of matches
      v2 = Ukk(sDiag,sCost, horz, d, c);
      if (v2>=0 && v2<lenA && A[v2] == B[v2-d] && (v2>res)) {
	res = v2; fromDir = horz; fromCost = c;
      }
      
      // Need to also check same cost for vert followed by a run of matches
      v3 = Ukk(sDiag,sCost, vert, d, c);
      if (v3>=0 && v3<lenA && A[v3] == B[v3-d] && (v3>res)) {
	res = v3; fromDir = vert; fromCost = c;
      }
      fromDiag = d;
	
      // Extend the diagonal if possible (only done for diagonal step)
      while (res>=0 && res<lenA && A[res] == B[res-d]) res++;
      break;
    case horz : 
      v1 = Ukk(sDiag,sCost, diag, d+1, c-a-b); 
      v2 = Ukk(sDiag,sCost, horz, d+1, c-b); 
      v3 = Ukk(sDiag,sCost, vert, d+1, c-a-b);
      switch (fromDir=MAXDIR3(v1,v2,v3)) {
      case diag : res = v1; ; fromCost = c-a-b; break;
      case horz : res = v2; ; fromCost = c-b; break;
      case vert : res = v3; ; fromCost = c-a-b; break;
      default :	assert(0);
      }
      fromDiag = d+1;
      
      break;
    case vert : 
      v1 = Ukk(sDiag,sCost, diag, d-1, c-a-b) + 1; 
      v2 = Ukk(sDiag,sCost, horz, d-1, c-a-b) + 1; 
      v3 = Ukk(sDiag,sCost, vert, d-1, c-b) + 1;
      switch (fromDir=MAXDIR3(v1,v2,v3)) {
      case diag : res = v1; fromCost = c-a-b; break;
      case horz : res = v2; fromCost = c-a-b; break;
      case vert : res = v3; fromCost = c-b; break;
      default :	assert(0);
      }
      fromDiag = d-1;
      
      break;
    default : assert(0);
    }

    if (res<0) res=BIG_NEGATIVE;

    data[d+offset][c%MODSIZE][matrix].val  = res; // Store the result in the 
    data[d+offset][c%MODSIZE][matrix].done = 1;   // data array
    data[d+offset][c%MODSIZE][matrix].cost = c;

    if (c<=checkpCost && c>checkpCost-(CHECKPSIZE-1)) {
      // This is on a checkpoint column (cost) so checkpoint it.
      storeCheckp(c, d, matrix, res, fromCost, fromDiag, fromDir);
    } else if (checkp[fromDiag+offset][fromCost%CHECKPSIZE][fromDir].cost>=0) {
      // Update checkpoint information if where step from had checkpoint info
      checkp[d+offset][c%CHECKPSIZE][matrix] =
	checkp[fromDiag+offset][fromCost%CHECKPSIZE][fromDir];
    }


    // Calculate all other entries for this cell (diag/cost).
    // This may do some unnessary work, but it simplifies the code.
    if (matrix!=diag) Ukk(sDiag, sCost, diag, d, c);
    if (matrix!=horz) Ukk(sDiag, sCost, horz, d, c);
    if (matrix!=vert) Ukk(sDiag, sCost, vert, d, c);
    
    
#ifdef DEBUG
    depth-=2;
#endif
    
    return res;			//Return the length obtainable on diag for cost
  }

  int doUkk(int sDiag, int sCost, direction sDir, int sDist,
	    int fDiag, int fCost, direction fDir)
  {
    int cost;
    int finalDist;

    if (fCost==sCost)
      return sDist;

    for (int i=0; i<arraySize;i++) // Initialize the main structure
      for (int j=0; j<MODSIZE; j++) 
	for (int k=0; k<NUMDIRS; k++)
	  data[i][j][k].done = 0;

    for (int i=0; i<arraySize; i++) // Initialize the checkpoint structure
      for (int j=0; j<CHECKPSIZE; j++)
	for (int k=0; k<NUMDIRS; k++)
	  checkp[i][j][k].cost = BIG_NEGATIVE;

    if (fDiag == sDiag)	     // Determine what cost calulation should start at
      cost = sCost + 1;
    else
      cost = sCost + abs(fDiag-sDiag);

    // Calculate what cost to checkpoint at.
    checkpCost = (fCost+cost+(CHECKPSIZE-1))/2;	

    data[sDiag+offset][sCost%MODSIZE][sDir].val  = sDist; //Fill in first point
    data[sDiag+offset][sCost%MODSIZE][sDir].cost = sCost;
    data[sDiag+offset][sCost%MODSIZE][sDir].done = 1;

    if (sCost<=checkpCost && sCost>checkpCost-(CHECKPSIZE-1))
      storeCheckp(sCost, sDiag, sDir, sDist, sCost, sDiag, sDir);
    
    do { // Main loop.

      // ----------------------------------------
      // Blank out the entries that will be calculated on this iteration
      // Corresponds to 2 diagonals from the point [diag][cost].
      for (int i=cost-MODSIZE,j=fDiag; i-sCost>=j-sDiag; i--,j++) {
	for (int k=0; k<NUMDIRS; k++)
	  data[j+offset][i%MODSIZE][k].done = 0;
      }
      for (int i=cost-MODSIZE,j=fDiag; i-sCost>=-(j-sDiag); i--,j--) {
	for (int k=0; k<NUMDIRS; k++)
	  data[j+offset][i%MODSIZE][k].done = 0;
      }
      // ----------------------------------------

      if (fDir==-1) {		// If don't know final dir calc them all
				// this is only done on the first time through
				// after that the exit direction is always
				// known (from the checkpoint data)
	int v2 = Ukk(sDiag, sCost, horz, fDiag, cost);
	int v3 = Ukk(sDiag, sCost, vert, fDiag, cost);
	int v1 = Ukk(sDiag, sCost, diag, fDiag, cost);
	finalDist = MAX3(v1,v2,v3);
	if (cost == fCost)
	  fDir = MAXDIR3(v1,v2,v3);
      } else
	finalDist = Ukk(sDiag, sCost, fDir, fDiag, cost);
      
      cost++;
    } while ((cost-1) != fCost);
    cost--;

    // Now do checkpoint recursion--------------------------------------------
    struct checkpElem cdata = checkp[fDiag+offset][fCost%CHECKPSIZE][fDir];
    if (cdata.cost<0) {
      fprintf(stderr,"INTERNAL ERROR: No data checkpointed\n");
      fprintf(stderr,"sCost=%d sDiag=%d sDir=%d  fCost=%d fDiag=%d fDir=%d"
	      " checkpCost=%d\n",sCost,sDiag,sDir,fCost,fDiag,fDir,checkpCost);
      exit(-1);
    }
    
    int r = doUkk(sDiag, sCost, sDir, sDist,
		  cdata.fromDiag, cdata.fromCost, cdata.fromDir);

    // Dislay one <ch1, ch2> with a possible diagonal extension
    char ch1,ch2;
    switch (cdata.entryDir) {
    case diag : ch1=A[r]; ch2=B[r-cdata.diag]; r++; break;
    case horz : ch1='-';  ch2=B[r-cdata.diag-1]; break;
    case vert : ch1=A[r]; ch2='-'; r++; break;
    default :	assert(0);
    }
    printf("<%c,%c> ",ch1,ch2);
    for(; r<cdata.dist; r++)
      printf("[%c,%c] ", A[r],B[r-cdata.diag]);

    // Now other half of recursion--------------------------------------------
    doUkk(cdata.diag, cdata.cost, cdata.entryDir, cdata.dist,
	  fDiag, fCost, fDir);
    

    return finalDist;
  }

public:
  // editCost() - Calculate edit distance and display the alignment between
  //              two strings A and B.
  int align(char strA[], char strB[], int fCost)
  {
    int cost;
    int finalDiag;
    int sDist;
    int maxEditdist;

    A = strA;
    B = strB;
    lenA = strlen(A);
    lenB = strlen(B);
    finalDiag = lenA-lenB;

    maxEditdist = (lenA>lenB) ? lenA : lenB; // Maximum possible edit distance.
    maxEditdist *= MODSIZE;	// This is an upper bound (I think)

    offset = (maxEditdist-finalDiag)/2; // Offset to allow for -ve diagonals
    
    // Now allocate storage. Largest possible edit distance by 3 columns needed
    data  = new struct ukkElem[maxEditdist+1][MODSIZE][NUMDIRS];
    arraySize = maxEditdist+1;

    checkp = new struct checkpElem[maxEditdist+1][CHECKPSIZE][NUMDIRS];

    
    // Calculate and display the initial matchings of A and B. (For entry 0,0
    // in the Ukkonen matrix.
    for (sDist=0; A[sDist] && A[sDist]==B[sDist]; sDist++) 
      printf("[%c,%c] ",A[sDist],B[sDist]);

    cost = doUkk(0, 0, diag, sDist, finalDiag, fCost, nodir);
    
    return cost;
  }
};

void msg(char *prog) {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit cost between two strings, and" << endl;
  cout << "displays an optimal alignment under linear gap costs.  This program uses a" << endl;
  cout << "modified version of Ukkonen's algorithm(1) and check-pointing(2) to recover the alignment." << endl;
  cout << "It has time complexity  O(n*log(d) + d*d) and space complexity O(d) (where d is the edit cost)" << endl;

  cout << endl;

  cout << "1:  E. Ukkonen, \"On Approximate String Matching\"," << endl;
  cout << "    Foundations of Computation Theory, 1983, 158, pp 487-495" << endl;

  cout << endl;

  cout << "2:  D. R. Powell, L. Allison and T. I. Dix," << endl;
  cout << "    \"A Versatile Divide and Conquer Technique for Optimal String Alignment\"," << endl;
  cout << "    Information Processing Letters, 1999, 70:3, pp 127-139" << endl;

  cout << endl;

  cout << "Usage: " << prog << endl;
  cout << endl << endl;
}


int main(int argc, char** argv)
{
  char A[MAXSTRING],B[MAXSTRING];
  int res;
  Ukkonen t;
  UkkonenCheckp t2;

  msg(argv[0]);

  printf("Match=%d Mis=%d a=%d b=%d\n",MatchCost,MismatchCost,a,b);

  Common::readStrings(A, B);

  res = t.editCost(A, B);
  t2.align(A,B,res);
  cout << endl << "Edit Cost = " << res << endl;
  return 0;
}






