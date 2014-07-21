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


// 
//
//
//   Date: 12/9/96
//   Author: David Powell

#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "common.h"

using namespace std;

#define a 3                     // insert(delete) cost = w(k) = a+b*k
#define b 1                     // Where k is number or inserts(deletes)
#define MatchCost    0          // Cost of a match (not used, must be 0)
#define MismatchCost 1          // Cost of a mismatch

				// Size of window on ukkonen matrix.
				// For a description of why this size is
				// necessary read window_size.note
#define MODSIZE      (MAX2((a+b)*2, MismatchCost)+1)
    

#define BIG_NEGATIVE -BIG_VAL

class Ukkonen
{
  struct ukkElem {
    int val;		     // Contents of this cell
    int column;              // Column of this cell (used for debugging)
    int done;		     // Flag for whether this cell has been calculated
  };
  
  enum direction {horz, vert, diag};
  
  char *A, *B;			// Two strings to be aligned
  int lenA,lenB;		// Length of the two strings
  struct ukkElem (*diagArray)[MODSIZE];
  struct ukkElem (*horzArray)[MODSIZE];
  struct ukkElem (*vertArray)[MODSIZE];
  
  int offset;			// Offset into data[] because -ve indices
				// needed.
  int arraySize;		// Number of rows allocated for each array

#if DEBUG
  int depth=-2;
#endif

private:
  int Ukk(int sDiag, int sCost, direction matrix, int d, int c)
  {
    int v1,v2,v3,res;
    struct ukkElem (*data)[MODSIZE];	// Pointer to the desired array

    switch (matrix) {
    case diag: data = diagArray; break;
    case horz: data = horzArray; break;
    case vert: data = vertArray; break;
    }
  
//    if (d==sDiag && c==sCost) return sDistance; // Starting condition

#if DEBUG
    for (v1=0;v1<depth+2;v1++)
      printf(" ");
    printf("D=%d C=%d DIR=%s\n",d,c,
	   (matrix==diag ? "Diag" : (matrix==horz ? "Horz" : "Vert")));
#endif
    
    
    if (abs(d-sDiag) > c-sCost)
      return BIG_NEGATIVE;	//Boundary condition

    if (data[d+offset][c%MODSIZE].done)	// Check if already calculated
      if (data[d+offset][c%MODSIZE].column==c)
        return data[d+offset][c%MODSIZE].val;
      else {
        fprintf(stderr,"INTERNAL ERROR: Use of column that is not correct");
        fprintf(stderr,"\nCost = %d,  Diag = %d,   (Cost in Matrix=%d)\n",
                c,d,data[d+offset][c%MODSIZE].column);
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
     case diag: 
      v1 = Ukk(sDiag,sCost, diag, d, c-MismatchCost)+1; 
      v2 = Ukk(sDiag,sCost, horz, d, c-MismatchCost)+1; 
      v3 = Ukk(sDiag,sCost, vert, d, c-MismatchCost)+1;
      res = MAX3(v1,v2,v3);


      v2 = Ukk(sDiag,sCost, horz, d, c);
      if (v2>=0 && v2<lenA && A[v2] == B[v2-d])
	res = MAX2(v2,res);
      
      v3 = Ukk(sDiag,sCost, vert, d, c);
      if (v3>=0 && v3<lenA && A[v3] == B[v3-d])
	res = MAX2(v3,res);
	
      // Extend the diagonal if possible (only done for diagonal step)
      while (res>=0 && res<lenA && A[res] == B[res-d]) res++;
      break;
    case horz: 
      v1 = Ukk(sDiag,sCost, diag, d+1, c-a-b); 
      v2 = Ukk(sDiag,sCost, horz, d+1, c-b); 
      v3 = Ukk(sDiag,sCost, vert, d+1, c-a-b);
      res = MAX3(v1,v2,v3);
      break;
    case vert: 
      v1 = Ukk(sDiag,sCost, diag, d-1, c-a-b) + 1; 
      v2 = Ukk(sDiag,sCost, horz, d-1, c-a-b) + 1; 
      v3 = Ukk(sDiag,sCost, vert, d-1, c-b) + 1;
      res = MAX3(v1,v2,v3);
      break;
    }

    if (res<0) res=BIG_NEGATIVE;

    data[d+offset][c%MODSIZE].val = res; // Store the result in the data array
    data[d+offset][c%MODSIZE].done= 1;
    data[d+offset][c%MODSIZE].column = c;

#ifdef DEBUG
    depth-=2;
#endif
    
    return res;			//Return the length obtainable on diag for cost
  }

  int doUkk(int sDiag, int sCost, int sDist, int fDiag)
  {
    int cost, dist;

    for (int i=0; i<arraySize;i++) { // Initialize the main structure
      for (int j=0; j<MODSIZE; j++) {
	diagArray[i][j].done = 0;
	horzArray[i][j].done = 0;
	vertArray[i][j].done = 0;
      }
    }

    if (fDiag == sDiag)	     // Determine what cost calulation should start at
      cost = sCost + 1;
    else
      cost = sCost + abs(fDiag-sDiag);
    
    diagArray[sDiag+offset][sCost%MODSIZE].val = sDist;
    diagArray[sDiag+offset][sCost%MODSIZE].done= 1;

    if (fDiag==sDiag) dist=sDist;
    else dist = -1;

    while (dist < lenA) { // Main loop.
 
      for (int i=cost-MODSIZE,j=fDiag; i>=abs(j); i--,j++) {
	diagArray[j+offset][i%MODSIZE].done = 0;
	horzArray[j+offset][i%MODSIZE].done = 0;
	vertArray[j+offset][i%MODSIZE].done = 0;
      }
      for (int i=cost-MODSIZE,j=fDiag; i>=abs(j); i--,j--) {
	diagArray[j+offset][i%MODSIZE].done = 0;
	horzArray[j+offset][i%MODSIZE].done = 0;
	vertArray[j+offset][i%MODSIZE].done = 0;
      }

      int v1 = Ukk(sDiag, sCost, diag, fDiag, cost);
      int v2 = Ukk(sDiag, sCost, horz, fDiag, cost);
      int v3 = Ukk(sDiag, sCost, vert, fDiag, cost);
      dist = MAX3(v1,v2,v3);
      
      cost++;
    };
    cost--;

    if (dist>lenA) {
      fprintf(stderr,"ERROR: Missed lenA.  Final dist=%d  lenA=%d\n",dist,lenA);
      exit(-1);
    }


    return cost;
  }

public:
  // editCost() - Calculate edit distance and display the alignment between
  //              two strings A and B.
  int editCost(char strA[], char strB[])
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
    diagArray = new struct ukkElem[maxEditdist+1][MODSIZE];
    horzArray = new struct ukkElem[maxEditdist+1][MODSIZE];
    vertArray = new struct ukkElem[maxEditdist+1][MODSIZE];
    arraySize = maxEditdist+1;
    
    // Calculate and display the initial matchings of A and B. (For entry 0,0
    // in the Ukkonen matrix.
    for (sDist=0; A[sDist] && A[sDist]==B[sDist]; sDist++) ;
//      printf("<%c,%c> ",A[sDist],B[sDist]);

    cost = doUkk(0, 0, sDist, finalDiag);
    
    return cost;
  }
};






