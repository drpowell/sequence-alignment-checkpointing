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


#ifndef __COMMON_H__
#define __COMMON_H__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#define MAXSTRING 20000          // Maximum size for reading in a string
#define BIG_VAL   1000000


#define MAX2(x,y) ((x)>(y) ? (x) : (y))
#define MAX3(x,y,z) ((x)>=(y) ? ((x)>=(z) ? (x) : (z)) : ((y)>=(z) ? (y) :(z)))

#define MIN3(x,y,z) ((x)<=(y) ? ((x)<=(z) ? (x) : (z)) : ((y)<=(z) ? (y) :(z)))



class Common {

public:

  static void readStrings(char *A, char *B) {

    cout << "Enter string A : ";
    if (!fgets(A, MAXSTRING, stdin)) {
      cerr << "Error reading input" << endl;
      exit(1);
    }

    cout << "Enter string B : ";
    if (!fgets(B, MAXSTRING, stdin)) {
      cerr << "Error reading input" << endl;
      exit(1);
    }

    A[strlen(A)-1] = 0;    // Lose the \n
    B[strlen(B)-1] = 0;
  }

};

#endif
