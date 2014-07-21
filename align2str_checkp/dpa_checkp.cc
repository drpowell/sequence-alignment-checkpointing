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


// DPA for 2 string alignment with checkpointing to find alignment.
// Time complexity=O(n*n)    Space complexity=O(n)
// Constant costs
//      match : 0
//      change,indel : 1
//   Date: 3/7/96
//   Author: David Powell
#include <ctype.h>
#include <iostream>
#include <stdio.h>

using namespace std;

#define PRINT

#define MAXSTRING 20000		// Maximum size for reading in a string

#define MIN3(x,y,z) ((x)<(y) ? ((x)<(z) ? (x) : (z)) : ((y)<(z) ? (y) : (z)))
#define MIN_INDEX3(x,y,z) ((x)<=(y) ? ((x)<=(z) ? 0 : 2) : ((y)<=(z) ? 1 : 2))

int alignment[2*MAXSTRING];
int alignPos = 0;

long loopCount;

struct crossingType {
  int splitPoint,exitPoint;
};

struct crossingType crossing[2][MAXSTRING+1];
int D[2][MAXSTRING+1];		// D is the table for alignment.
				// D should be sized dynamically. But I am too
				// lazy.

// doDpa -
int doDpa(char A[], char B[], int a0, int b0, int a1, int b1)
{
  int splitRow, n=a1-a0, m=b1-b0;

#ifdef DEBUG
  printf("Doing A[%d..%d] and B[%d..%d]\n",a0,a1,b0,b1);
#endif

  if (a0==a1 || b0==b1) {	// Is it a trivial case? 
				// (just inserts or deletes)
    if (a0==a1) {		// Is it inserts?
      for (int i=b0;i<b1;i++)
	alignment[alignPos++] = 1; // 1 == Insert
    } else {			// Must be deletes.
      for (int i=a0;i<a1;i++)
	alignment[alignPos++] = 2; // 1 == Delete
    }
    return 0;
  }
  
  splitRow = n/2;
  D[0][0] = 0;			// Initialize D
  for (int j=0;j<=m;j++)
    D[0][j] = j;
  
  for (int k=0;k<=m;k++) {	// Initialise crossing info. (only used if 
    crossing[0][k].splitPoint = k; // splitRow == 0).
    crossing[0][k].exitPoint = k-1;
  }
  
  // Calculate the D array for the edit distance
  for (int i=1;i<=n;i++) {
    D[i%2][0] = i;
    crossing[i%2][0].splitPoint = 0; crossing[i%2][0].exitPoint = 0;
    for (int j=1;j<=m;j++) {
      int matchCost, insertCost, deleteCost;

      loopCount++;
      
      matchCost  = D[(i-1)%2][j-1] + (A[i-1+(a0)]==B[j-1+(b0)] ? 0 : 1);
      insertCost = D[i%2][j-1] + 1;
      deleteCost = D[(i-1)%2][j] + 1;
      
      switch (MIN_INDEX3(matchCost, insertCost, deleteCost)) {
      case 0:			// Match or mismatch
	D[i%2][j] = matchCost;
	crossing[i%2][j] = crossing[(i-1)%2][j-1];
	if (i==splitRow+1)
	  crossing[i%2][j].exitPoint = j;
	break;
      case 1:			// Insert
	D[i%2][j] = insertCost;
	crossing[i%2][j] = crossing[i%2][j-1];
	// Note: no exit dir recorded if it was an insert
	//       and crossing info comes from same row.
	break;
      case 2:			// Delete
	D[i%2][j] = deleteCost;
	crossing[i%2][j] = crossing[(i-1)%2][j];
	if (i==splitRow+1)
	  crossing[i%2][j].exitPoint = j;
	break;
      }
    } // end for j

    // Set up check point if it is half way
    if (i==splitRow) {
      for (int k=0;k<=m;k++) {
	crossing[i%2][k].splitPoint = k;
	crossing[i%2][k].exitPoint = k-1; // Assume exit is Insert until
				          // determined on next pass
      }
    }
    
  } // end for i

  int editDistance = D[n%2][m];	// Save the actual edit distance.
  int splitColumn = crossing[n%2][m].splitPoint; // Where to finish top half
  int startPoint  = crossing[n%2][m].exitPoint;  // Where to start bottom half
  
				// Recurse for top half
  doDpa(A, B, a0, b0, a0+splitRow,  b0+splitColumn);
  
				// Now store alignment info. It was either a
				// delete or a match (mismatch), cause it had
				// to move to a new row.
  alignment[alignPos++] = (splitColumn==startPoint) ? 2 : 0;

				// Recurse for bottom half
  doDpa(A, B, a0+splitRow+1, b0+startPoint, a1, b1);
    
//  cout << "Edit distance = " << D[n%2][m] << endl;
//  cout << "Split point = " << splitRow << "," << splitColumn << endl;
//  cout << "Starting point for bottom half = " << startPoint << endl;
    
  return editDistance;		// Return edit distance
}

void msg() {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit distance between two strings, and" << endl;
  cout << "displays an optimal alignment.  This program uses a basic DPA with" << endl;
  cout << "check-pointing(1) to recover the alignment, and has time complexity O(n*n)," << endl;
  cout << "and space complexity O(n)" << endl;
  cout << endl;
  cout << "1:  D. R. Powell, L. Allison and T. I. Dix," << endl;
  cout << "    \"A Versatile Divide and Conquer Technique for Optimal String Alignment\"," << endl;
  cout << "    Information Processing Letters, 1999, 70:3, pp 127-139" << endl;

  cout << endl << endl;
}

int main(int argc, char *argv[])
{
  char A[MAXSTRING],B[MAXSTRING];
  int res;

  msg();

  cout << "Enter string A : ";
  cin >> A;
  cout << "Enter string B : ";
  cin >> B;

  // Ensure Strings are all uppercase
  for (int i=0;A[i];i++)
    A[i] = toupper(A[i]);
  for (int i=0;B[i];i++)
    B[i] = toupper(B[i]);

  loopCount = 0;
  res = doDpa(A, B, 0, 0, strlen(A), strlen(B));

#ifdef PRINT
  for (int i=0,j=0,pos=0;pos<alignPos;pos++) {
    int ch1,ch2;
    switch (alignment[pos]) {
    case 0:			// Match/mismatch
      ch1 = A[i++];
      ch2 = B[j++];
      break;
    case 1:			// Insertion
      ch1 = '-';
      ch2 = B[j++];
      break;
    case 2:			// Deletion
      ch1 = A[i++];
      ch2 = '-';
      break;
    }
    printf("<%c,%c> ",ch1,ch2);
  }
#endif
  cout << endl;
  cout << "Edit distance = " << res << endl;
  cout << "Loop Counter = " << loopCount << endl;

  return 0;
}







