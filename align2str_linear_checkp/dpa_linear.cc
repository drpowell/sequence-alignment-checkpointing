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

// Does edit distance and alignment of two string using standard DPA and linear
// insert/delete costs.

#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"

using namespace std;


#define _a 3			// insert(delete) cost = w(k) = a+b*k
#define _b 1			// Where k is number or inserts(deletes)
#define _MatchCost    0		// Cost of a match
#define _MismatchCost 1		// Cost of a mismatch

int a = _a;
int b = _b;
int MatchCost = _MatchCost;
int MismatchCost = _MismatchCost;

class DPAlinear
{
private:

  
  struct dpaElem {
    int horz;
    int vert;
    int diag;
  };

  typedef struct dpaElem* dpaElemPtr;

  enum direction {horz, vert, diag};
  
  char *A,*B;			// The two strings being compared

  int lenA,lenB;
  
  struct dpaElem **data;	// Where the 2D DPA array will reside.
				// Uses an array of pointers to each row.

private:

  direction minDirection(struct dpaElem d, int i, int j)
  {
    direction res;
    if (i==0) d.vert = BIG_VAL;
    if (j==0) d.horz = BIG_VAL;
    
    if (d.horz <= d.vert)
      if (d.horz <= d.diag)
	res = horz;
      else
	res = diag;
    else if (d.vert <= d.diag)
      res = vert;
    else
      res = diag;
    return res;
  }

  void dispAlign()
  {
    int i,j;
    char alignment[2*(lenA+lenB)];
    int pos=0;
    
    i = lenA;
    j = lenB;
    direction dir = minDirection(data[i][j],i,j);

    while (i!=0 || j!=0) {
      char ch1, ch2;
      struct dpaElem d;
      
      ch1 = ch2 = '-';
      switch (dir) {
      case horz:
	ch2 = B[j-1];
	j--;
	d.diag = data[i][j].diag + a + b;
	d.horz = data[i][j].horz +     b;
	d.vert = data[i][j].vert + a + b;
	break;
      case vert:
	ch1 = A[i-1];
	i--;
	d.diag = data[i][j].diag + a + b;
	d.horz = data[i][j].horz + a + b;
	d.vert = data[i][j].vert +     b;
	break;
      case diag:
	ch1 = A[i-1];
	ch2 = B[j-1];
	i--;
	j--;
	d.diag = data[i][j].diag;
	d.horz = data[i][j].horz;
	d.vert = data[i][j].vert;
	break;
      }

      dir = minDirection(d,1,1);
      alignment[pos++] = ch1;
      alignment[pos++] = ch2;
    }

    for (i=pos-1; i>0; i-=2)
      printf("<%c,%c> ",alignment[i-1],alignment[i]);
  }
  
  int doDPA()
  {
    int i,j;

    for (i=1;i<=lenA;i++) {
      data[i][0].horz = BIG_VAL;
      data[i][0].vert = a+b*i;
      data[i][0].diag = BIG_VAL;
    }
    for (j=1;j<=lenB;j++) {
      data[0][j].horz = a+b*j;
      data[0][j].vert = BIG_VAL;
      data[0][j].diag = BIG_VAL;
    }
    data[0][0].diag = 0;
    data[0][0].vert = BIG_VAL;
    data[0][0].horz = BIG_VAL;

    for (i=1; i<=lenA; i++) {
//      printf("\n%d  ",i);
      for (j=1; j<=lenB; j++) {

	data[i][j].horz = MIN3(
	  data[i][j-1].horz + b,
	  data[i][j-1].vert + a + b,
	  data[i][j-1].diag + a + b);

	data[i][j].vert = MIN3(
	  data[i-1][j].horz + a + b,
	  data[i-1][j].vert + b,
	  data[i-1][j].diag + a + b);
	
	int diagCost = (A[i-1]==B[j-1]) ? MatchCost : MismatchCost;
	
	data[i][j].diag = MIN3(
	  data[i-1][j-1].horz + diagCost,
	  data[i-1][j-1].vert + diagCost,
	  data[i-1][j-1].diag + diagCost);
	
//        printf("H=%d V=%d D=%d    ",data[i][j].horz,data[i][j].vert,data[i][j].diag);
      }
    }

    dispAlign();

    return MIN3(data[lenA][lenB].horz,
		data[lenA][lenB].vert,
		data[lenA][lenB].diag);
  }
  
public:
  int doAlign(char strA[], char strB[])
  {
    struct dpaElem *tmp;
    int i;
    
    A = strA;
    B = strB;
    lenA = strlen(A);
    lenB = strlen(B);

    data = new dpaElemPtr[lenA+1];
    tmp = new struct dpaElem[(lenA+1) * (lenB+1)];
    for (i=0;i<=lenA;i++)
      data[i] = &tmp[i * (lenB+1)];

    return doDPA(); 
  }
  
};

void msg(char *prog) {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit cost between two strings, and" << endl;
  cout << "displays an optimal alignment under linear gap costs.  This program uses a" << endl;
  cout << "basic DPA and has time and space complexity of O(n*n)" << endl;
  cout << endl;
  cout << "Usage: " << prog << " [matchCost mismatchCost a b]" << endl;
  cout << "  where cost for gap of length k = a + b*k" << endl;
  cout << endl << endl;
}

int main(int argc, char *argv[])
{
  char A[MAXSTRING],B[MAXSTRING];
  int cost;
  DPAlinear t;

  msg(argv[0]);
  
  if (argc==5) {
    MatchCost = atoi(argv[1]);
    MismatchCost = atoi(argv[2]);
    a = atoi(argv[3]);
    b = atoi(argv[4]);
  }

  printf("Match=%d Mis=%d a=%d b=%d\n",MatchCost,MismatchCost,a,b);

  Common::readStrings(A, B);

  cost = t.doAlign(A,B);

  cout << "Edit cost = " << cost << endl;
  

  return 0;
}

