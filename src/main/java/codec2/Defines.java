/*
 * Copyright (C) 1990-2013 David Rowe, VK5DGRe
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
 * Global definitions
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public class Defines {

    public static final float TWO_PI = (float) (2.0 * Math.PI);
    //
    public static final int FFT_SIZE = 512;          // size of FFT used
    public static final int FS = 8000;               // sample rate in Hz
    public static final int N = 80;                  // number of samples per 10 ms frame
    public static final int M = 320;                 // pitch analysis frame size
    public static final int WO_E_BITS = 8;
    public static final int WO_BITS = 7;
    public static final int E_BITS = 5;
    public static final int WO_LEVELS = (1 << WO_BITS);
    public static final int E_LEVELS = (1 << E_BITS);
    public static final int P_MIN = 20;              // minimum pitch
    public static final int P_MAX = 160;             // maximum pitch
    public static final float WO_MAX = TWO_PI / P_MIN;     // MAX as in the bigger number
    public static final float WO_MIN = TWO_PI / P_MAX;     // MIN as in the smaller number
    public static final float STEP = (WO_MAX - WO_MIN) / WO_LEVELS;
    public static final int LPC_ORD = 10;
    public static final int NW = 279;                // analysis window size
}
