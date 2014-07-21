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


// file: ukk_checkp2
// This file uses Ukkonen's alogorithm to align 2 strings in O(d) space by
// using checkpoints on columns of the ukkonen matrix. Then recursion is
// performed on the two parts of the matrix (before and after the checkpoint).
// This is similar to Hirshberg's technique for the DPA.

#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "ukk_noalign.h"

using namespace std;

#define MAXSTRING 20000          // Maximum size for reading in a string
#define BIG_NEGATIVE -10        // Well maybe not that big :)

#define INSERT     -1           // Direction of insert. (decrement diagonal)
#define MISMATCH   0            // Same diagonal
#define DELETE     1            // Increment diagonal

#define MAX_DIR3(x,y,z) ((x)>(y) ? ((x)>=(z) ? INSERT : DELETE) : ((y)>=(z) ? MISMATCH : DELETE)) // MAX_DIR3 - Determine whether to insert, delete, mismatch.


#define MAX2(x,y) ((x)>(y) ? (x) : (y))
#define MIN2(x,y) ((x)<(y) ? (x) : (y))
#define MAX3(x,y,z) ((x)>(y) ? ((x)>=(z) ? (x) : (z)) : ((y)>=(z) ? (y) :(z)))

class Ukkonen_align
{
private:
  
  struct checkpointT {		//Structure to store the checkpoint information
    int diagonal;		// Diagonal that this cell comes from.
    int distance;		// The cell contents (distance along A)
    int entryDir;		// Direction for entry into this cell.
  };
  
  char *A,*B;			// The two strings being compared
  
  int (*data)[2];		// Array for ukkonen's algorithm
  struct checkpointT (*cpData)[2]; // checkpoint data

  int offset;			// offset into the arrary (for -ve diagonals)

  int numInsert, numMismatch, numDelete, numMatch; // some statistics variables
  
public:
  int innerLoop,outerLoop;

private:
  // do_Ukk - Performs Ukkonens alogorithm from a give diagonal, cost to
  // a finishing diagonal, cost. The ukkonen matrix is checkpointed at a cost
  // half way between the starting and finishing cost. Then the do_Ukk recurses
  // to determine the alignment.
  int do_Ukk(int sDiag, int sCost, int sDist, int fDiag, int fCost)
  {
    int checkpoint,checkpointed;

    if (sCost==fCost)		// Trivial case.
      return sDist;

    checkpoint = (fCost-sCost+1)/2 + sCost; // Cost to checkpoint at.
    checkpointed = 0;

    data[sDiag+offset][sCost%2] = sDist; // Starting point of the alogorithm
     
    for (int cost=sCost+1; cost<=fCost; cost++) {
      int i;

      // Now loop over the diagonals.
      // The MAX2/MIN2 calculation determines the necessary starting/finishing
      // point so no unnessary work is done.
      for (i = MAX2(sDiag-(cost-sCost), fDiag-(fCost-cost));
	   i <= MIN2(sDiag+(cost-sCost), fDiag+(fCost-cost));
	   i++) {
	int v1,v2,v3,res;
	
	v1 = v2 = v3 = BIG_NEGATIVE;
	if (i+1 <= sDiag+(cost-1-sCost))
	  v1 = data[i+1+offset][(cost-1)%2];
	if (i >= sDiag-(cost-1-sCost) && i <= sDiag+(cost-1-sCost))
	  v2 = data[i+offset][(cost-1)%2]+1;
	if (i-1 >= sDiag-(cost-1-sCost))
	  v3 = data[i-1+offset][(cost-1)%2]+1;
	
	res = MAX3(v1,v2,v3);

	while (A[res] && A[res] == B[res-i]) {// Extend the diagonal
	  res++;
	  innerLoop++;
	}
	
	data[i+offset][cost%2] = res;

	int dir = MAX_DIR3(v1,v2,v3); // Was it an insert, delete or mismatch?
	
	if (checkpointed) {	// If checkpointed already, move CP info along
	  cpData[i+offset][cost%2] = cpData[i+offset-dir][(cost-1)%2];
	} else {             // If not, keep CP data up to date for when we do.
	  cpData[i+offset][cost%2].diagonal = i;
	  cpData[i+offset][cost%2].distance = res;
	  cpData[i+offset][cost%2].entryDir = dir;
	}
	outerLoop++;
      } //for
      
      if (cost == checkpoint)
	checkpointed = 1;
    }

    // Now the must determine the exact matrix cell to split at.
    int splitCost = checkpoint;
    int splitDiag = cpData[fDiag+offset][fCost%2].diagonal;
    int splitDist = cpData[fDiag+offset][fCost%2].distance;
    int splitDir  = cpData[fDiag+offset][fCost%2].entryDir;
    int fDist = data[fDiag+offset][fCost%2];

    int r = do_Ukk(sDiag, sCost, sDist, splitDiag-splitDir, splitCost-1);

    dispAlignment(splitDir, splitDiag-splitDir, r, splitDist);

    do_Ukk(splitDiag, splitCost, splitDist, fDiag, fCost);
    
    return fDist;        // Return distance reached for finishing cost/diagonal
  }
  
  // dispAlignment() - Display a single non-match (ie. mismatch, insertion or
  // deletion) followed by a number of matches.
  void dispAlignment(int dir, int diag, int dist, int newDist)
  {
    char ch1='-',ch2='-';
    switch (dir) {
    case INSERT:
      ch2 = B[dist-diag];
      numInsert++;
      break;
    case MISMATCH:
      ch1 = A[dist];
      ch2 = B[dist-diag];
      numMismatch++;
      break;
    case DELETE:
      ch1 = A[dist];
      numDelete++;
      break;
    }
    printf("<%c,%c> ",ch1,ch2);
    
    if (dir!=INSERT)
      dist++;
    for (int i=0;i<newDist-dist;i++, numMatch++)
     printf("<%c,%c> ",A[i+dist],B[i+dist-diag-dir]);
  }    

  
public:
  // doAlign - is the entry point for the string alignment. It is necessary
  // to also supply the edit distance for the two strings.
  void doAlign(char strA[], char strB[], int editDist)
  {
    int lenA,lenB;
    int sDist,maxCost;

    A = strA;
    B = strB;
    lenA = strlen(A);
    lenB = strlen(B);

    numMatch = numMismatch = numInsert = numDelete = 0;
    innerLoop = outerLoop = 0;

    
    maxCost = MAX2(lenA,lenB);
    offset = (maxCost-(lenA-lenB))/2;
    
    data = new int[maxCost+1][2]; // Create room for Ukkonen matrix
    cpData = new struct checkpointT[maxCost+1][2]; // Make the CP data array.

    // Calculate and display the initial matchings of A and B. (For entry 0,0
    // in the Ukkonen matrix.
    for (sDist=0; A[sDist] && A[sDist]==B[sDist]; sDist++,numMatch++)
      printf("<%c,%c> ",A[sDist],B[sDist]);
    
    do_Ukk(0,0,sDist,lenA-lenB,editDist);

    cout << endl;
  }
  
};

void msg() {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit distance between two strings, and displays" << endl;
  cout << "an optimal alignment.  This program uses Ukkonen's algorithm(1) with " << endl;
  cout << "check-pointing(2) to recover the alignment.  The average time complexity" << endl;
  cout << "is O(n*log(d) + d*d), and space complexity is O(d)   (where d is the edit distance)" << endl;

  cout << endl;

  cout << "1:  E. Ukkonen, \"On Approximate String Matching\"," << endl;
  cout << "    Foundations of Computation Theory, 1983, 158, pp 487-495" << endl;

  cout << endl;

  cout << "2:  D. R. Powell, L. Allison and T. I. Dix," << endl;
  cout << "    \"A Versatile Divide and Conquer Technique for Optimal String Alignment\"," << endl;
  cout << "    Information Processing Letters, 1999, 70:3, pp 127-139" << endl;

  cout << endl << endl;

}

int main(int argc, char *argv[])
{
  char A[MAXSTRING],B[MAXSTRING];
  int cost;
  Ukkonen_align t;
  Ukkonen t2;

  msg();

  cout << "Enter string A : ";
  cin >> A;
  cout << "Enter string B : ";
  cin >> B;

  
  cost = t2.editCost(A,B);	// First determined the edit distance
  
  cout << endl << "LenA="<<strlen(A)<<"   lenB="<<strlen(B);
  cout << endl << "Base: Inner = " << t2.innerLoop
       << "   Outer = " << t2.outerLoop << endl;

  cout << "Edit distance = " << cost << endl;
 
  t.doAlign(A,B,cost);		// Now determine the alignment.
  
  cout << endl << "Align: Inner = "<< t.innerLoop
       << "   Outer = " << t.outerLoop << endl;

  return 0;
}

