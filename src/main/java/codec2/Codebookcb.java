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

public final class Codebookcb {

    public static Codebook[] lsp_cb;

    public Codebookcb() {
        lsp_cb = new Codebook[10];

        /*
         * 36 bits total, and 132 vector elements
         */
        lsp_cb[0] = new Codebook(1, 4, 16, codes0);
        lsp_cb[1] = new Codebook(1, 4, 16, codes1);
        lsp_cb[2] = new Codebook(1, 4, 16, codes2);
        lsp_cb[3] = new Codebook(1, 4, 16, codes3);
        lsp_cb[4] = new Codebook(1, 4, 16, codes4);
        lsp_cb[5] = new Codebook(1, 4, 16, codes5);
        lsp_cb[6] = new Codebook(1, 4, 16, codes6);
        lsp_cb[7] = new Codebook(1, 3, 8, codes7);
        lsp_cb[8] = new Codebook(1, 3, 8, codes8);
        lsp_cb[9] = new Codebook(1, 2, 4, codes9);
    }
    private final static float codes0[] = {
        225.0F,
        250.0F,
        275.0F,
        300.0F,
        325.0F,
        350.0F,
        375.0F,
        400.0F,
        425.0F,
        450.0F,
        475.0F,
        500.0F,
        525.0F,
        550.0F,
        575.0F,
        600.0F
    };
    private final static float codes1[] = {
        325.0F,
        350.0F,
        375.0F,
        400.0F,
        425.0F,
        450.0F,
        475.0F,
        500.0F,
        525.0F,
        550.0F,
        575.0F,
        600.0F,
        625.0F,
        650.0F,
        675.0F,
        700.0F
    };
    private final static float codes2[] = {
        500.0F,
        550.0F,
        600.0F,
        650.0F,
        700.0F,
        750.0F,
        800.0F,
        850.0F,
        900.0F,
        950.0F,
        1000.0F,
        1050.0F,
        1100.0F,
        1150.0F,
        1200.0F,
        1250.0F
    };
    private final static float codes3[] = {
        700.0F,
        800.0F,
        900.0F,
        1000.0F,
        1100.0F,
        1200.0F,
        1300.0F,
        1400.0F,
        1500.0F,
        1600.0F,
        1700.0F,
        1800.0F,
        1900.0F,
        2000.0F,
        2100.0F,
        2200.0F
    };
    private final static float codes4[] = {
        950.0F,
        1050.0F,
        1150.0F,
        1250.0F,
        1350.0F,
        1450.0F,
        1550.0F,
        1650.0F,
        1750.0F,
        1850.0F,
        1950.0F,
        2050.0F,
        2150.0F,
        2250.0F,
        2350.0F,
        2450.0F
    };
    private final static float codes5[] = {
        1100.0F,
        1200.0F,
        1300.0F,
        1400.0F,
        1500.0F,
        1600.0F,
        1700.0F,
        1800.0F,
        1900.0F,
        2000.0F,
        2100.0F,
        2200.0F,
        2300.0F,
        2400.0F,
        2500.0F,
        2600.0F
    };
    private final static float codes6[] = {
        1500.0F,
        1600.0F,
        1700.0F,
        1800.0F,
        1900.0F,
        2000.0F,
        2100.0F,
        2200.0F,
        2300.0F,
        2400.0F,
        2500.0F,
        2600.0F,
        2700.0F,
        2800.0F,
        2900.0F,
        3000.0F
    };
    private final static float codes7[] = {
        2300.0F,
        2400.0F,
        2500.0F,
        2600.0F,
        2700.0F,
        2800.0F,
        2900.0F,
        3000.0F
    };
    private final static float codes8[] = {
        2500.0F,
        2600.0F,
        2700.0F,
        2800.0F,
        2900.0F,
        3000.0F,
        3100.0F,
        3200.0F
    };
    private final static float codes9[] = {
        2900.0F,
        3100.0F,
        3300.0F,
        3500.0F
    };
}