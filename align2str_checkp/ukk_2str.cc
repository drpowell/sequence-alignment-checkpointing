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


// Ukkonnen's algorithm for 2 string edit distance (alignment not determined)
// costs:   match = 0,   change = indel = 1
//   Date: 15/7/96
//   Author: David Powell
#include <ctype.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "ukk_noalign.h"

using namespace std;


#define MAXSTRING 20000		// Maximum size for reading in a string

void msg() {
  cout << "Copyright (C) David Powell <david@drp.id.au>" << endl;
  cout << "  This program comes with ABSOLUTELY NO WARRANTY; and is provided" << endl;
  cout << "  under the GNU Public License v2, for details see file COPYRIGHT" << endl << endl;

  cout << "This program calculates the edit distance between two strings, _but_ does not" << endl;
  cout << "recover an alignment.  This program uses Ukkonen's algorithm(1) and has" << endl;
  cout << "average time complexity of O(d*d + n), and space complexity O(d)   (where d is the edit distance)" << endl;

  cout << endl;

  cout << "1:  E. Ukkonen, \"On Approximate String Matching\"," << endl;
  cout << "    Foundations of Computation Theory, 1983, 158, pp 487-495" << endl;

  cout << endl << endl;

}

int main(int argc, char *argv[])
{
  char A[MAXSTRING],B[MAXSTRING];
  int res;
  Ukkonen t;

  msg();

  cout << "Enter string A : ";
  cin >> A;
  cout << "Enter string B : ";
  cin >> B;

  res = t.editCost(A,B);
  cout << endl << "Edit distance = " << res << endl;
  return 0;
}






