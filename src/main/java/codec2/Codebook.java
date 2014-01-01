/*
 * Copyright (C) 1990-2013 David Rowe, VK5DGR
 *
 * All Rights Reserved
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 2.1, as published
 * by the Free Software Foundation. This program is distributed in the hope that
 * it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 */
package codec2;

/**
 * Class that describes each codebook. More like a C structure actually.
 *
 * <p>Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Codebook {

    public int k;         // dimension of the vector
    public int log2m;     // log base 2 of size m
    public int m;         // number of vector elements
    public float[] cb;    // the actual codebook array

    public Codebook(int val1, int val2, int val3, float[] val4) {
        this.k = val1;
        this.log2m = val2;
        this.m = val3;
        this.cb = val4;
    }
}
