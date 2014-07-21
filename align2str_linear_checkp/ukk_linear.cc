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
#include "ukk_linear.h"

using namespace std;

void msg(char *prog) {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit cost between two strings under linear gap costs, _but_" << endl;
  cout << "does determine an alignment.  This program uses a modified version of Ukkonen's algorithm(1)" << endl;
  cout << "and has time complexity O(d*d + n) and space complexity O(d) (where d is the edit cost)" << endl;

  cout << endl;

  cout << "1:  E. Ukkonen, \"On Approximate String Matching\"," << endl;
  cout << "    Foundations of Computation Theory, 1983, 158, pp 487-495" << endl;

  cout << endl;

  cout << "Usage: " << prog << endl;
  cout << endl << endl;
}

int main(int argc, char **argv)
{
  char A[MAXSTRING],B[MAXSTRING];
  int res;
  Ukkonen t;

  msg(argv[0]);

  printf("Match=%d Mis=%d a=%d b=%d\n",MatchCost,MismatchCost,a,b);

  Common::readStrings(A, B);

  res = t.editCost(A,B);
  cout << "Edit cost = " << res << endl;
  return 0;
}







