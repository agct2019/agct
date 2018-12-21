package test;

/*
 * Copyright (c) 2005 by L. Paul Chew.
 *
 * Permission is hereby granted, without written agreement and without
 * license or royalty fees, to use, copy, modify, and distribute this
 * software and its documentation for any purpose, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

/**
 * Points in Euclidean space, implemented as double[].
 * <p/>
 * Includes simple geometric operations. Uses matrices; a matrix is represented
 * as an array of Pnts. Uses simplices; a simplex is represented as an array of
 * Pnts.
 *
 * @author Paul Chew
 *         <p/>
 *         Created July 2005. Derived from an earlier, messier version.
 */
public class Pnt {
    /*
     * // The point's coordinates
     */
    public double[] coordinates;
    /*
     * // The gene index in the selected genes
     */
    private int geneIndex;

    /*
     * Constructor.
     *
     * @param coords the coordinates
     */
    public Pnt(double[] coords, int name) {
        // Copying is done here to ensure that Pnt's coords cannot be altered.
        coordinates = new double[coords.length];
        System.arraycopy(coords, 0, coordinates, 0, coords.length);
        geneIndex = name;
    }

    public Pnt(Pnt copy) {
        this(copy.coordinates, copy.geneIndex);
    }

    /**
     * d次元の点を宣言する。 coordinates はdoubleの配列
     */
    public Pnt(int d) {
        coordinates = new double[d];
    }

    /**
     * Constructor.
     *
     * @param coordA
     * @param coordB
     * @param coordC
     */
    public Pnt(double coordA, double coordB, double coordC, int name) {
        this(new double[]{coordA, coordB, coordC}, name);
    }

    public int getGeneIndex() {
        return geneIndex;
    }

    /**
     * Create a String for this Pnt.
     *
     * @return a String representation of this Pnt.
     */
    public String toString() {
        if (coordinates.length == 0)
            return "()";
        String result = "Pnt[" + geneIndex + "](" + coordinates[0];
        for (int i = 1; i < coordinates.length; i++)
            result = result + "," + coordinates[i];
        result = result + ")";
        return result;
    }

    /**
     * Equality.
     *
     * @param other the other Object to compare to
     * @return true iff the Pnts have the same coordinates
     */
    public boolean equals(Object other) {
        if (!(other instanceof Pnt))
            return false;
        Pnt p = (Pnt) other;
        if (this.coordinates.length != p.coordinates.length)
            return false;
        for (int i = 0; i < this.coordinates.length; i++)
            if (this.coordinates[i] != p.coordinates[i])
                return false;
        return true;
    }

    /**
     * HashCode.
     *
     * @return the hashCode for this Pnt
     */
    public int hashCode() {
        int hash = 0;
        for (int i = 0; i < this.coordinates.length; i++) {
            long bits = Double.doubleToLongBits(this.coordinates[i]);
            hash = (31 * hash) ^ (int) (bits ^ (bits >> 32));
        }
        return hash;
    }

    /* Pnts as vectors */

    /**
     * @return this Pnt's dimension.
     */
    public int dimension() {
        return coordinates.length;
    }

    /**
     * Check that dimensions match.
     *
     * @param p the Pnt to check (against this Pnt)
     * @return the dimension of the Pnts
     * @throws IllegalArgumentException if dimension fail to match
     */
    public int dimCheck(Pnt p) {
        int len = this.coordinates.length;
        if (len != p.coordinates.length)
            throw new IllegalArgumentException("Dimension mismatch");
        return len;
    }

    /**
     * Create a new Pnt by adding additional coordinates to this Pnt.
     *
     * @param coords the new coordinates (added on the right end)
     * @return a new Pnt with the additional coordinates
     */
    public Pnt extend(double[] coords) {
        double[] result = new double[coordinates.length + coords.length];
        System.arraycopy(coordinates, 0, result, 0, coordinates.length);
        System.arraycopy(coords, 0, result, coordinates.length, coords.length);
        return new Pnt(result, geneIndex);
    }

    /**
     * Dot product.
     *
     * @param p the other Pnt
     * @return dot product of this Pnt and p
     */
    public double dot(Pnt p) {
        int len = dimCheck(p);
        double sum = 0;
        for (int i = 0; i < len; i++)
            sum += this.coordinates[i] * p.coordinates[i];
        return sum;
    }

    /**
     * Magnitude (as a vector).
     *
     * @return the Euclidean length of this vector
     */
    public double magnitude() {
        return Math.sqrt(this.dot(this));
    }

    /* Pnts as matrices */

    /**
     * Create a String for a matrix.
     *
     * @param matrix the matrix (an array of Pnts)
     * @return a String represenation of the matrix
     */
    public static String toString(Pnt[] matrix) {
        StringBuffer buf = new StringBuffer("{");
        for (int i = 0; i < matrix.length; i++)
            buf.append(" " + matrix[i]);
        buf.append(" }");
        return buf.toString();
    }

    /**
     * Compute the determinant of a matrix (array of Pnts). This is not an
     * efficient implementation, but should be adequate for low dimension.
     *
     * @param matrix the matrix as an array of Pnts
     * @return the determinnant of the input matrix
     * @throws IllegalArgumentException if dimensions are wrong
     */
    public static double determinant(Pnt[] matrix) {
        if (matrix.length != matrix[0].dimension())
            throw new IllegalArgumentException("Matrix is not square");
        boolean[] columns = new boolean[matrix.length];
        for (int i = 0; i < matrix.length; i++)
            columns[i] = true;
        try {
            return determinant(matrix, 0, columns);
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IllegalArgumentException("Matrix is wrong shape");
        }
    }

    /**
     * Compute the determinant of a submatrix specified by starting row and by
     * "active" columns.
     *
     * @param matrix  the matrix as an array of Pnts
     * @param row     the starting row
     * @param columns a boolean array indicating the "active" columns
     * @return the determinant of the specified submatrix
     * @throws ArrayIndexOutOfBoundsException if dimensions are wrong
     */
    private static double determinant(Pnt[] matrix, int row, boolean[] columns) {
        if (row == matrix.length)
            return 1;
        double sum = 0;
        int sign = 1;
        for (int col = 0; col < columns.length; col++) {
            if (!columns[col])
                continue;
            columns[col] = false;
            sum += sign * matrix[row].coordinates[col]
                    * determinant(matrix, row + 1, columns);
            columns[col] = true;
            sign = -sign;
        }
        return sum;
    }

    /**
     * Compute generalized cross-product of the rows of a matrix. The result is
     * a Pnt perpendicular (as a vector) to each row of the matrix. This is not
     * an efficient implementation, but should be adequate for low dimension.
     *
     * @param matrix the matrix of Pnts (one less row than the Pnt dimension)
     * @return a Pnt perpendicular to each row Pnt
     * @throws IllegalArgumentException if matrix is wrong shape
     */
    public static Pnt cross(Pnt[] matrix) {
        int len = matrix.length + 1;
        if (len != matrix[0].dimension())
            throw new IllegalArgumentException("Dimension mismatch");
        boolean[] columns = new boolean[len];
        for (int i = 0; i < len; i++)
            columns[i] = true;
        double[] result = new double[len];
        int sign = 1;
        try {
            for (int i = 0; i < len; i++) {
                columns[i] = false;
                result[i] = sign * determinant(matrix, 0, columns);
                columns[i] = true;
                sign = -sign;
            }
        } catch (ArrayIndexOutOfBoundsException e) {
            throw new IllegalArgumentException("Matrix is wrong shape");
        }
        return new Pnt(result, -1);
    }

    /* Pnts as simplices */

    /**
     * Determine the signed content (i.e., area or volume, etc.) of a simplex.
     *
     * @param simplex the simplex (as an array of Pnts)
     * @return the signed content of the simplex
     */
    public static double content(Pnt[] simplex) {
        Pnt[] matrix = new Pnt[simplex.length];
        for (int i = 0; i < matrix.length; i++)
            matrix[i] = simplex[i].extend(new double[]{1});
        int fact = 1;
        for (int i = 1; i < matrix.length; i++)
            fact = fact * i;
        return determinant(matrix) / fact;
    }

    /**
     * Relation between this Pnt and a simplex (represented as an array of
     * Pnts). Result is an array of signs, one for each vertex of the simplex,
     * indicating the relation between the vertex, the vertex's opposite facet,
     * and this Pnt.
     * <p/>
     * <pre>
     *   -1 means Pnt is on same side of facet
     *    0 means Pnt is on the facet
     *   +1 means Pnt is on opposite side of facet
     * </pre>
     *
     * @param simplex an array of Pnts representing a simplex
     * @return an array of signs showing relation between this Pnt and the
     *         simplex
     * @throws IllegalArgumentExcpetion if the simplex is degenerate
     */
    public int[] relation(Pnt[] simplex) {
    /*
     * In 2D, we compute the cross of this matrix: 1 1 1 1 p0 a0 b0 c0 p1 a1
	 * b1 c1 where (a, b, c) is the simplex and p is this Pnt. The result is
	 * a vector in which the first coordinate is the signed area (all signed
	 * areas are off by the same constant factor) of the simplex and the
	 * remaining coordinates are the *negated* signed areas for the
	 * simplices in which p is substituted for each of the vertices.
	 * Analogous results occur in higher dimensions.
	 */
        int dim = simplex.length - 1;
        if (this.dimension() != dim)
            throw new IllegalArgumentException("Dimension mismatch");

	/* Create and load the matrix */
        Pnt[] matrix = new Pnt[dim + 1];
	/* First row */
        double[] coords = new double[dim + 2];
        for (int j = 0; j < coords.length; j++)
            coords[j] = 1;
        matrix[0] = new Pnt(coords, -1);
	/* Other rows */
        for (int i = 0; i < dim; i++) {
            coords[0] = this.coordinates[i];
            for (int j = 0; j < simplex.length; j++)
                coords[j + 1] = simplex[j].coordinates[i];
            matrix[i + 1] = new Pnt(coords, -1);
        }

	/* Compute and analyze the vector of areas/volumes/contents */
        Pnt vector = cross(matrix);
        double content = vector.coordinates[0];
        int[] result = new int[dim + 1];
        for (int i = 0; i < result.length; i++) {
            double value = vector.coordinates[i + 1];
            if (Math.abs(value) <= 1.0e-6 * Math.abs(content))
                result[i] = 0;
            else if (value < 0)
                result[i] = -1;
            else
                result[i] = 1;
        }
        if (content < 0) {
            for (int i = 0; i < result.length; i++)
                result[i] = -result[i];
        }
        if (content == 0) {
            for (int i = 0; i < result.length; i++)
                result[i] = Math.abs(result[i]);
        }
        return result;
    }

    /**
     * Test if this Pnt is outside of simplex.
     *
     * @param simplex the simplex (an array of Pnts)
     * @return the simplex Pnt that "witnesses" outsideness (or null if not
     *         outside)
     */
    public Pnt isOutside(Pnt[] simplex) {
        int[] result = this.relation(simplex);
        for (int i = 0; i < result.length; i++) {
            if (result[i] > 0)
                return simplex[i];
        }
        return null;
    }

    /**
     * Test relation between this Pnt and circumcircle of a simplex.
     *
     * @param simplex the simplex (as an array of Pnts)
     * @return -1, 0, or +1 for inside, on, or outside of circumcircle
     */
    public int vsCircumcircle(Pnt[] simplex) {
        Pnt[] matrix = new Pnt[simplex.length + 1];
        for (int i = 0; i < simplex.length; i++)
            matrix[i] = simplex[i].extend(new double[]{1,
                    simplex[i].dot(simplex[i])});
        matrix[simplex.length] = this
                .extend(new double[]{1, this.dot(this)});
        double d = determinant(matrix);
        int result = (d < 0) ? -1 : ((d > 0) ? +1 : 0);
        if (content(simplex) < 0)
            result = -result;
        return result;
    }

    public Pnt subtract(Pnt x) {
        if (coordinates.length != x.coordinates.length)
            Matrix.perror("Pnt.class :: coordinates mistmatch");

        Pnt a = new Pnt(x.coordinates.length);
        for (int i = 0; i < x.coordinates.length; i++) {
            a.coordinates[i] = coordinates[i] - x.coordinates[i];
        }
        return a;
    }

    public boolean all_null() {
        int i;
        for (i = 0; i < coordinates.length; i++)
            if (coordinates[i] != 0.0)
                return false;
        return true;
    }

}