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


// Does edit distance and alignment of 2 strings with linear insert/delete 
// costs.  Alignment is recovered by using checkpointing (similar to Hirschberg)
// So Time complexity O(n*n)    Space complexity O(n)
// (Not counting local variables at each recursion step (neither does 
//  Hirschberg in his analysis).


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

#define MINDIR3(h,v,d) ((h)<=(v) ? \
                        ((h)<=(d) ? horz : diag) : \
                        ((v)<=(d) ? vert : diag))

int a = _a;
int b = _b;
int MatchCost = _MatchCost;
int MismatchCost = _MismatchCost;


class DPAlinear
{
  private:

    // enumerated type for the possible directions
    enum direction {any=-1, horz=0, vert=1, diag=2};

    // a 'state' corresponds to the incoming direction in
    // the dpa matrix
    struct stateType {
      int cost;

      int fromCol;
      direction fromDir;
    };

    // Each cell of the dpa matrix has 3 states, one for
    // each possible previous opertation
    struct dpaElem {
      struct stateType d[3];
    };

    char *A,*B;			// The two strings being compared

    int lenA,lenB;

    struct dpaElem *data[2];	// Where the 2D DPA array will reside.

    // Store the final alignment here - in reverse!
    char *alignment;
    int alignPos;

    int debugPrint;

  private:
 
    // cellCalc(state, from, addCosts[]) - critical routine
    // that fills in 1 'state' of a cell in the DPA matrix.
    // 'from' is the determining cell which contains 3
    // states.  addCosts[] is a array of 3 integer costs
    // for the cost of transition from a state to this one.
    // The three costs in this matrix are only different
    // when the transition is from a horz->horz or
    // vert->vert state, that is, continuing a gap.
    // This routine also copies along the 'from' info, so
    // the crossing column, and state, of the check-point
    // row can be determined.
    void cellCalc(struct stateType &state, const struct dpaElem &from, int addCosts[3]) {
      int c[3];
      for (int i=0; i<3; i++)
        c[i] = from.d[i].cost + addCosts[i];

      int i = MINDIR3(c[horz], c[vert], c[diag]);

      state.cost = c[i];
      state.fromCol = from.d[i].fromCol;
      state.fromDir = from.d[i].fromDir;
    }

    // The DPA with check-pointing recursion to recover an
    // alignment.  Note: on first call eDir==any because the
    // final state of the alignment is not constrained.  All
    // subsequent calls will have eDir set to a proper
    // value.  sDir==diag on the first call.
    int doDPA(
        int start_i, int start_j,   direction sDir,
        int finish_i, int finish_j, direction eDir)
    {
      int checkRow;   // Index of row to check-point
      int i,j;

      checkRow  = (finish_i - start_i + 1)/2 + start_i;

      // Initialise starting state (based on sDir)
      for (int dir=0; dir<3; dir++)
        data[start_i%2][start_j].d[dir].cost = (dir == sDir ? 0 : BIG_VAL);

      // Iterate over the DPA matrix
      for (i=start_i; i<=finish_i; i++) {
        for (j=start_j; j<=finish_j; j++) {

          int addCosts[3];

          if (i==start_i && j==start_j)
            continue;    // Start state already initialised.

          // Calculate the horizontal state (delete a char from strA)
          if (j>start_j) {
            addCosts[horz] = b;                      // Continue a horizontal gap
            addCosts[vert] = addCosts[diag] = a+b;   // Start a horizontal gap
            cellCalc(data[i%2][j].d[horz], data[i%2][j-1], addCosts);
          } else
            data[i%2][j].d[horz].cost = BIG_VAL;

          // Calculate the vertical state (delete a char from strB)
          if (i>start_i) {
            addCosts[vert] = b;                    // Continue a vertical gap
            addCosts[horz] = addCosts[diag] = a+b; // Start a vertical gap
            cellCalc(data[i%2][j].d[vert], data[(i-1)%2][j], addCosts);
          } else
            data[i%2][j].d[vert].cost = BIG_VAL;

          // Calculate the diagonal state (match/mismatch)
          if (i>start_i && j>start_j) {
            addCosts[horz] = (A[i-1]==B[j-1]) ? MatchCost : MismatchCost;
            addCosts[vert] = addCosts[diag] = addCosts[horz];
            cellCalc(data[i%2][j].d[diag], data[(i-1)%2][j-1], addCosts);
          } else
            data[i%2][j].d[diag].cost = BIG_VAL;
    
          if (debugPrint) printf("\n(%d,%d) H=%d V=%d D=%d", i,j,data[i%2][j].d[horz].cost,data[i%2][j].d[vert].cost,data[i%2][j].d[diag].cost);

          // If this is the check-point row, store the 'from' info
          if (i==checkRow)
            for (int dir=0; dir<3; dir++) {
              data[i%2][j].d[dir].fromCol = j;
              data[i%2][j].d[dir].fromDir = (direction)dir;
            }
        }
      }


      if (eDir < 0)  // Determine the end state if our caller did not set it.
        eDir = MINDIR3(data[finish_i%2][finish_j].d[horz].cost,
                       data[finish_i%2][finish_j].d[vert].cost,
                       data[finish_i%2][finish_j].d[diag].cost);

      int editDist = data[finish_i%2][finish_j].d[eDir].cost;  // Final edit cost

      if (finish_i - start_i > 1) {
        // More than 2 rows in matrix.  Recurse to find alignment
        int split_j = data[finish_i%2][finish_j].d[eDir].fromCol;
        direction split_dir = data[finish_i%2][finish_j].d[eDir].fromDir;

        doDPA(checkRow, split_j, split_dir, finish_i, finish_j, eDir);
        doDPA(start_i, start_j, sDir, checkRow, split_j, split_dir);
      } else {
        // Determine alignment directly from 'data'.  Only 1 or 2 rows
        i = finish_i;
        j = finish_j;
        direction dir = eDir;
        int all_horz = (i == start_i);  // all_horz - true when only 1 row in this step
        while (i != start_i || j != start_j) {
          int val;
          char c1,c2;
          switch (dir) {
          case vert:
            c1 = A[i-1]; c2 = '-';
            i--;
            dir = horz;
            all_horz = 1;
            break; 
          case diag:
            c1 = A[i-1]; c2 = B[j-1];
            i--; j--;
            dir = horz;
            all_horz = 1;
            break; 
          case horz:
            c1 = '-'; c2 = B[j-1];
            val = data[i%2][j].d[horz].cost;
            j--;

            if (all_horz || val == data[i%2][j].d[horz].cost + b)
              dir = horz;
            else if (val == data[i%2][j].d[vert].cost + a + b)
              dir = vert;
            else
              dir = diag;
            break;
          case any:
            printf("BAD logic!");  exit(911);
          }
          alignment[alignPos++] = c1;
          alignment[alignPos++] = c2;
        }
      }

      return editDist;
    }

    int doCheckpDPA ()  
    {
      int i,res;

      //debugPrint=1;
      res = doDPA(0, 0, diag, lenA, lenB, any);

      for (i=alignPos-2; i>=0; i-=2)
        printf("<%c,%c> ", alignment[i], alignment[i+1]);
  
      return res;
    }
  
public:
  int doAlign(char strA[], char strB[])
  {
    struct dpaElem *tmp;
    
    A = strA;
    B = strB;
    lenA = strlen(A);
    lenB = strlen(B);

    tmp = new struct dpaElem[2 * (lenB+1)];
    data[0] = &tmp[0];
    data[1] = &tmp[1*(lenB+1)];

    alignment = new char[(lenA+lenB)*2];
    alignPos = 0;
    
    return doCheckpDPA();
  }
  
};

void msg(char *prog) {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit cost between two strings, and" << endl;
  cout << "displays an optimal alignment under linear gap costs.  This program uses a" << endl;
  cout << "basic DPA and check-pointing(1) to recover the alignment.  It has time complexity" << endl;
  cout << "of O(n*n), and space complexity of O(n)" << endl;
  cout << endl;
  cout << "1:  D. R. Powell, L. Allison and T. I. Dix," << endl;
  cout << "    \"A Versatile Divide and Conquer Technique for Optimal String Alignment\"," << endl;
  cout << "    Information Processing Letters, 1999, 70:3, pp 127-139" << endl;

  cout << endl << endl;

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

