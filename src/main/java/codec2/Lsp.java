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
 * Line Spectral Pairs (LSP) and Linear Predictive Coding (LPC) algorithms
 *
 * This class contains methods for LPC to LSP conversion and LSP to LPC
 * conversion.
 *
 * The LSP coefficients are not in radians format but in the x domain of the
 * unit circle.
 *
 * <p>Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Lsp {

    private final float[] Table;    
    
    public Lsp() {
        Table = new float[(Defines.LPC_ORD / 2) + 1];              // 0..5
        Table[0] = 1.0F;                                           // first value always 1.0
    }
    
    /**
     * This method evaluates a series of chebyshev polynomials
     *
     * @param coef a float array coefficients of the polynomial to be evaluated
     * @param x a float point where polynomial is to be evaluated
     * @return a float representing the polynomial evaluation
     */
    private float cheb_poly_eval(float[] coef, float x) {
        int m = (Defines.LPC_ORD / 2);        // 5
        int i;

        /*
         * Evaluate chebyshev series formulation using iterative approach
         */
        Table[1] = x;

        for (i = 2; i <= m; i++) {          // 2..5
            Table[i] = 2.0F * x * Table[i - 1] - Table[i - 2];
        }

        /*
         * Evaluate polynomial and return value
         */

        float sum = 0.0F;

        for (i = 0; i <= m; i++) {          // 0..5
            sum += (coef[m - i] * Table[i]);
        }

        return sum;
    }

    /**
     * This method converts LPC coefficients to LSP coefficients.
     *
     * @param a float array representing the lpc coefficients
     * @param freq a float array representing the LSP frequencies in radians
     * @param nb an int representing the number of sub-intervals
     * @param delta a float representing the grid spacing interval
     * @return an int representing number of roots found
     */
    public int lpc_to_lsp(float[] a, float[] freq, int nb, float delta) {
        float[] Q = new float[Defines.LPC_ORD + 1];
        float[] P = new float[Defines.LPC_ORD + 1];
        float psuml, psumr, psumm, temp_xr, xl, xr, xm, temp_psumr;
        int i, j, k;
        float[] px;                	// ptrs of respective P'(z) & Q'(z)
        int px_i = 0;
        float[] qx;
        int qx_i = 0;
        float[] p;
        int p_i = 0;
        float[] q;
        int q_i = 0;
        float[] pt;                	// ptr used for cheb_poly_eval()
        //int pt_i = 0;                 // whether P' or Q'
        int roots = 0;                  // number of roots found
        boolean flag = true;
        int m = Defines.LPC_ORD / 2;    // order of P'(z) & Q'(z) polynimials

        /*
         * determine P'(z)'s and Q'(z)'s coefficients where
         * P'(z) = P(z)/(1 + z^(-1)) and Q'(z) = Q(z)/(1-z^(-1))
         */

        px = p = P;                     // initialise ptrs
        qx = q = Q;

        px[px_i++] = 1.0F;
        qx[qx_i++] = 1.0F;

        for (i = 1; i <= m; i++) {
            px[px_i++] = a[i] + a[Defines.LPC_ORD + 1 - i] - p[p_i++];
            qx[qx_i++] = a[i] - a[Defines.LPC_ORD + 1 - i] + q[q_i++];
        }

        px = P;
        qx = Q;
        px_i = qx_i = 0;

        for (i = 0; i < m; i++) {
            px[px_i] = 2.0F * px[px_i];
            qx[qx_i] = 2.0F * qx[qx_i];
            px_i++;
            qx_i++;
        }

        /*
         * Search for a zero in P'(z) polynomial first and then alternate to Q'(z).
         * Keep alternating between the two polynomials as each zero is found
         */

        xr = xm = 0.0F;           	// initialise xr and xm to zero
        xl = 1.0F;               	// start at point xl = 1

        for (j = 0; j < Defines.LPC_ORD; j++) {
            if ((j % 2) != 0) {         // determines whether P' or Q' is eval.
                pt = Q;
            } else {
                pt = P;
            }

            psuml = cheb_poly_eval(pt, xl);          // evals poly. at xl
            flag = true;

            while (flag && (xr >= -1.0)) {
                xr = xl - delta;                    // interval spacing
                psumr = cheb_poly_eval(pt, xr);      // poly(xl-delta_x)
                temp_psumr = psumr;
                temp_xr = xr;

                /*
                 * if no sign change increment xr and re-evaluate
                 * poly(xr). Repeat til sign change.  if a sign change has
                 * occurred the interval is bisected and then checked again
                 * for a sign change which determines in which interval the
                 * zero lies in.  If there is no sign change between poly(xm)
                 * and poly(xl) set interval between xm and xr else set
                 * interval between xl and xr and repeat till root is located
                 * within the specified limits
                 */

                if ((psumr * psuml) < 0.0) {
                    roots++;

                    psumm = psuml;

                    for (k = 0; k <= nb; k++) {
                        xm = (xl + xr) / 2.0F;              // bisect the interval
                        psumm = cheb_poly_eval(pt, xm);

                        if (psumm * psuml > 0.0) {
                            psuml = psumm;
                            xl = xm;
                        } else {
                            psumr = psumm;
                            xr = xm;
                        }
                    }

                    // once zero is found, reset initial interval to xr

                    freq[j] = xm;
                    xl = xm;
                    flag = false;                 // reset flag for next search
                } else {
                    psuml = temp_psumr;
                    xl = temp_xr;
                }
            }
        }

        // convert from x domain to radians

        for (i = 0; i < Defines.LPC_ORD; i++) {
            freq[i] =  (float) Math.acos(freq[i]);
        }

        return roots;
    }

    /**
     * This method converts LSP coefficients to LPC coefficients.
     *
     * @param lsp array of LSP frequencies in radians
     * @param ak array of LPC coefficients
     */
    public void lsp_to_lpc(float[] lsp, float[] ak) {
        int i, j;
        float xout1, xout2, xin1, xin2;
        float[] pw;
        int pw_i = 0;       // simulate pointers in Java
        float[] n1;
        int n1_i = 0;
        float[] n2;
        int n2_i = 0;
        float[] n3;
        int n3_i = 0;
        float[] n4;
        int n4_i = 0;

        int m = Defines.LPC_ORD / 2;
        float[] freq = new float[Defines.LPC_ORD];
        float[] Wp = new float[(Defines.LPC_ORD * 4) + 2];  // java inits to 0.0

        // convert from radians to the x=cos(w) domain

        for (i = 0; i < Defines.LPC_ORD; i++) {
            freq[i] =  (float) Math.cos(lsp[i]);
        }

        // Set pointers up

        pw = n1 = n2 = n3 = n4 = Wp;
        xin1 = 1.0F;
        xin2 = 1.0F;

        /*
         * reconstruct P(z) and Q(z) by cascading second order polynomials
         * in form 1 - 2xz(-1) +z(-2), where x is the LSP coefficient
         */

        for (j = 0; j <= Defines.LPC_ORD; j++) {
            for (i = 0; i < m; i++) {
                n1_i = pw_i + (i * 4);
                n2_i = n1_i + 1;
                n3_i = n2_i + 1;
                n4_i = n3_i + 1;

                xout1 = xin1 - 2.0F * (freq[2 * i]) * n1[n1_i] + n2[n2_i];
                xout2 = xin2 - 2.0F * (freq[2 * i + 1]) * n3[n3_i] + n4[n4_i];

                n2[n2_i] = n1[n1_i];
                n4[n4_i] = n3[n3_i];
                n1[n1_i] = xin1;
                n3[n3_i] = xin2;

                xin1 = xout1;
                xin2 = xout2;
            }

            xout1 = xin1 + n4[n4_i + 1];
            xout2 = xin2 - n4[n4_i + 2];
            ak[j] = (xout1 + xout2) * 0.5F;

            n4[n4_i + 1] = xin1;
            n4[n4_i + 2] = xin2;

            xin1 = 0.0F;
            xin2 = 0.0F;
        }
    }
}