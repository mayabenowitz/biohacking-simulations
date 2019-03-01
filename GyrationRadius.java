/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import ffx.potential.bonded.Atom;
import java.util.logging.Logger;

/**
 *
 * @author Jbenovitz
 */
// WARNING! The system's center of mass must be set to zero for each timestep.
// Numbers in [] correspond to line numbers in the script.
// [25] Initializes a 3x3 matrix calculated by taking the outer product of the postion vector with itself.
// [26] Initializes a column position vector.
// [27] Initializes a row position vector, the transpose.
// [28] k: Iterates over the 1st atom in atomArray to the last atom.
// [29] i: Iterates over the columns of a 3x3 matrix and the column position vector.
// [30] j: Iterates over the rows of a 3x3 matrix and the row position vector.
// [31] For every atom in atomArray, the (x,y,z) coordinates are stored at the ith location in the column position vector.
// [32] For every atom in atomArray, the (x,y,z) coordinates are stored at the jth location in the row position vector.
// [33] For every (x,y,z) coordinate stored in the column and position vector, every (i,j) pair is multiplied. This is the outer product.
public class GyrationRadius {

    private static final Logger logger = Logger.getLogger(GyrationRadius.class.getName());
    // This is a property that GyrationRadius can have.
    private Atom atomArray[];

    // This creates an object of type "GyrationRadius"
    public GyrationRadius(Atom allAtoms[]) {
        atomArray = allAtoms;
        double[][] Outer_Product = new double[3][3];
        double[] Position_Vector = new double[3];
        double[] Position_Vector_Transpose = new double[3];
        for (int k = 0; k < atomArray.length; k++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Position_Vector[i] = atomArray[k].getXYZ(null)[i];
                    Position_Vector_Transpose[j] = atomArray[k].getXYZ(null)[j];
                    Outer_Product[i][j] = Position_Vector[i] * Position_Vector_Transpose[j];
                    logger.info(String.format("Outer_Product[%d][%d] = %f", i, j, Outer_Product[i][j]));
                }
            }
        }
    }

    private void generateGyrationMatrix() {
    }
}
