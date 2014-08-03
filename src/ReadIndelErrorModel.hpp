/*
 * Dindel
 * http://www.sanger.ac.uk/resources/software/dindel/
 *
 * Copyright 2010, Kees Albers
 * Released under the GPL Version 3.
 */
/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * ReadIndelErrorModel.hpp
 *
 *  Created on: Oct 12, 2009
 *      Author: caa
 */

#ifndef READINDELERRORMODEL_HPP_
#define READINDELERRORMODEL_HPP_

#include <vector>

class ReadIndelErrorModel {
public:
    ReadIndelErrorModel() {
        double hp[] = { 2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3 };
        baselineProbs = std::vector<double>(10, 0.0);
        for (int x = 0; x < 10; x++) {
            baselineProbs[x] = hp[x];
        }
    }

    double getViterbiHPError(int hpLen) const {
        int len = hpLen;
        if (len < 1) {
            len = 1;
        }
        double pbe;
        if (len <= 10) {
            pbe = baselineProbs[len - 1];
        } else {
            pbe = baselineProbs[9] + 4.3e-4 * double(len - 10);
        }
        pbe *= double(hpLen);
        if (pbe > 0.99) {
            pbe = 0.99;
        }
        return pbe;
    }

private:
    std::vector<double> baselineProbs;
};

#endif /* READINDELERRORMODEL_HPP_ */
