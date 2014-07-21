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


// The most basic DPA for 2 string alignment with costs :
//      match : 0
//      change,indel : 1
//   Date: 2/7/96
//   Author: David Powell
#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <malloc.h>

using namespace std;

#define MAXSTRING 20000		// Maximum size for reading in a string

#define MIN3(x,y,z) ((x)<(y) ? ((x)<(z) ? (x) : (z)) : ((y)<(z) ? (y) : (z)))

// dispArray - display the DPA array
void dispArray(int D[], char A[], char B[])
{
  int n,m;
  n = strlen(A);
  m = strlen(B);
  cout << endl << "     ";
  for (int j=0;j<=m;j++)
    cout << (j>0 ? B[j-1] : ' ') << "   ";
  for (int i=0;i<=n;i++) {
    cout << endl << (i>0 ? A[i-1] : ' ') << "   ";
    for (int j=0;j<=m;j++)
      printf("%2d  ",D[i*(m+1)+j]);
  }
  cout << endl;
}

#define P(i,j) ((i)*(m+1)+(j))	// P is a macro to access 1 dimensional D as
				// if it were 2 dimensional

// dispAlignment - Backtraces on the D array to find the actual alignment.
void dispAlignment(int D[], char A[], char B[])
{
  int n,m,i,j,pos;
  char chars[2*MAXSTRING],ch1,ch2;
  n = strlen(A);
  m = strlen(B);
  i=n;
  j=m;
  pos=0;
  // Work out alignment in reverse order, then display in correct order
  // Choice of alignment out of all the optimal alignments makes choice in this
  // order: Match, Insert, Delete   (on the reverse alignment)
  while (i!=0 || j!=0) {
    ch1 = ch2 = '-';
    if (i==0) {			// Left column, must be insert
      ch2 = B[j-1];
      j--;
    } else if (j==0) {		// Top row, must be delete
      ch1 = A[i-1];
      i--;
    } else {
      int matchCost, insertCost, deleteCost;
      matchCost  = D[P(i-1,j-1)] + (A[i-1]==B[j-1] ? 0 : 1);
      insertCost = D[P(i,j-1)] + 1;
      deleteCost = D[P(i-1,j)] + 1;
      if (matchCost <= insertCost) {
	if (matchCost <= deleteCost) { 
	  ch1 = A[i-1];		// A match or mismatch
	  ch2 = B[j-1];
	  i--;
	  j--;
	} else {
	  ch1 = A[i-1];		// A Delete
	  i--;
	}
      } else if (insertCost <= deleteCost) {
	ch2 = B[j-1];		// An Insert
	j--;
      } else {
	ch1 = A[i-1];		// A Delete
	i--;
      }
    }
    chars[pos++] = ch1;
    chars[pos++] = ch2;
  }

  // Display the alignment
  for (pos-=2;pos>=0;pos-=2)
    printf("<%c,%c> ",chars[pos],chars[pos+1]);
  cout << endl;
}

long loopCount;

// doDpa - Does standard DPA on 2 strings A and B.
int doDpa(char A[], char B[])
{
  int res,m,n;
  n = strlen(A);
  m = strlen(B);

  // Ensure Strings are all uppercase
  for (int i=0;i<n;i++)
    A[i] = toupper(A[i]);
  for (int i=0;i<m;i++)
    B[i] = toupper(B[i]);
  {
    int *Data;
    Data = (int *)malloc(sizeof(int)*(n+1)*(m+1));
    if (!Data) {fprintf(stderr,"Unable to malloc memory\n"); exit(-1); }

#define D(x,y) Data[(x)*(m+1)+y]
    
    D(0,0) = 0;		// Initialize D
    for (int i=1;i<=n;i++)
      D(i,0) = i;
    for (int j=1;j<=m;j++)
      D(0,j) = j;
    
    // Calculate the D array for the edit distance
    for (int i=1;i<=n;i++)
      for (int j=1;j<=m;j++) {
        loopCount++;
	D(i,j) = MIN3(
	  D(i-1,j-1) + (A[i-1]==B[j-1] ? 0 : 1), // Match or change
	  D(i,j-1) + 1,	// Insert
	  D(i-1,j) + 1);	// Delete
      }
    
    // Display the DPA array
//    dispArray((int*)D,A,B);

    // Display the alignment
    dispAlignment(Data,A,B);
    
    res = D(n,m);		// Store edit distance
    free(Data);
  }
  return res;
}

void msg() {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit distance between two strings, and" << endl;
  cout << "displays an optimal alignment.  This program uses a basic DPA and has" << endl;
  cout << "time and space complexity of O(n*n)" << endl;
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

  loopCount = 0;
  res = doDpa(A,B);
  cout << endl;
  cout << "Edit distance = " << res << endl;
  cout << "Loop Counter = " << loopCount << endl;
  return 0;
}






