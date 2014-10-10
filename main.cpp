/* Copyright (C) 2014  Ward Poelmans

This file is part of DOCI-Exact.

Doci-Exact is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Hubbard-GPU is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Hubbard-GPU.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <boost/timer.hpp>
#include <getopt.h>

#include "Permutation.h"
#include "Molecule.h"
#include "Hamiltonian.h"

/**
 * @mainpage
 * This is an exact DOCI solver by means of a lanczos solver. We build the
 * Hamiltonian and store it as a sparse matrix (in Column Compressed Format).
 * The Permutation object generates the next basis state on the fly.
 * The Molecule object holds the atomic/molecular integrals to use.
 *
 * @author Ward Poelmans <wpoely86@gmail.com>
 * @version   0.1
 * @date      2014
 * @copyright GNU Public License v3
 */

int main(int argc, char **argv)
{
    using std::cout;
    using std::endl;

    cout.precision(10);

    std::string integralsfile = "mo-integrals.h5";

    struct option long_options[] =
    {
        {"integrals",  required_argument, 0, 'i'},
        {"help",  no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int i,j;

    while( (j = getopt_long (argc, argv, "hi:", long_options, &i)) != -1)
        switch(j)
        {
            case 'h':
            case '?':
                cout << "Usage: " << argv[0] << " [OPTIONS]\n"
                    "\n"
                    "    -i, --integrals=integrals-file  Set the input integrals file\n"
                    "    -h, --help                      Display this help\n"
                    "\n";
                return 0;
                break;
            case 'i':
                integralsfile = optarg;
                break;
        }

    cout << "Reading: " << integralsfile << endl;

    Molecule mol(integralsfile);




    return 0;
}

/* vim: set ts=8 sw=4 tw=0 expandtab :*/
