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
 * Non Linear Pitch (NLP) estimation
 *
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Nlp {

    private static final float PREV_BIN = (4000.0F / (float) Math.PI) * (Defines.FFT_SIZE * Nlp.DEC) / Defines.FS;
    private static final float CNLP = 0.3F;          // post processor constant
    private static final float COEFF = 0.95F;        // notch filter parameter
    private static final int NLP_NTAP = 48;          // Decimation LPF order
    private static final int PMAX_M = 320;           // maximum NLP analysis window size
    private static final int DEC = 5;                // decimation factor
    private static final int MIN_BIN = Defines.FFT_SIZE * DEC / Defines.P_MAX;   // 16
    //
    private final FFT fftEncode;
    private final Complex[] fw;                      // DFT of squared signal (input)
    //
    private final float[] mem_fir;                   // decimation FIR filter memory
    private final float[] sq;                        // squared speech samples
    private final float[] coswintab;                 // optimization
    private float mem_x;                             // memory for notch filter
    private float mem_y;
    private float pitch;

    /*
     * 48 tap 600Hz low pass FIR filter coefficients
     * (symmetric around 24 points)
     */
    private final float[] nlp_fir = {
        -1.0818124e-03F,
        -1.1008344e-03F,
        -9.2768838e-04F,
        -4.2289438e-04F,
        5.5034190e-04F,
        2.0029849e-03F,
        3.7058509e-03F,
        5.1449415e-03F,
        5.5924666e-03F,
        4.3036754e-03F,
        8.0284511e-04F,
        -4.8204610e-03F,
        -1.1705810e-02F,
        -1.8199275e-02F,
        -2.2065282e-02F,
        -2.0920610e-02F,
        -1.2808831e-02F,
        3.2204775e-03F,
        2.6683811e-02F,
        5.5520624e-02F,
        8.6305944e-02F,
        1.1480192e-01F,
        1.3674206e-01F,
        1.4867556e-01F,
        1.4867556e-01F,
        1.3674206e-01F,
        1.1480192e-01F,
        8.6305944e-02F,
        5.5520624e-02F,
        2.6683811e-02F,
        3.2204775e-03F,
        -1.2808831e-02F,
        -2.0920610e-02F,
        -2.2065282e-02F,
        -1.8199275e-02F,
        -1.1705810e-02F,
        -4.8204610e-03F,
        8.0284511e-04F,
        4.3036754e-03F,
        5.5924666e-03F,
        5.1449415e-03F,
        3.7058509e-03F,
        2.0029849e-03F,
        5.5034190e-04F,
        -4.2289438e-04F,
        -9.2768838e-04F,
        -1.1008344e-03F,
        -1.0818124e-03F
    };

    public Nlp(FFT fft) {
        fftEncode = fft;

        this.fw = new Complex[Defines.FFT_SIZE];
        this.mem_fir = new float[NLP_NTAP];
        this.sq = new float[PMAX_M];                // java initializes to zero
        this.mem_x = 0.0F;
        this.mem_y = 0.0F;
        this.pitch = 0.0F;
        this.coswintab = new float[Defines.M / Nlp.DEC];

        /*
         * Pre-compute decimation window
         */
        for (int i = 0; i < Defines.M / Nlp.DEC; i++) {
            coswintab[i] = (float) (0.5 - 0.5 * Math.cos(2.0 * Math.PI * (double) i / (double) (Defines.M / Nlp.DEC - 1)));
        }
    }

    public float getPitch() {
        return this.pitch;
    }

    /**
     * *
     * Called by analyze_one_frame() in Codec2
     *
     * <p>
     * Determines the pitch in samples using a Non Linear Pitch (NLP) algorithm.
     *
     * <p>
     * Returns the fundamental in Hz. Note that the actual pitch estimate is for
     * the center of the M sample Sn[] vector, not the current N sample input
     * vector. This is (I think) a delay of 2.5 frames with N=80 samples. You
     * should align further analysis using this pitch estimate to be centered on
     * the middle of Sn[].
     *
     * <p>
     * In the presence of background noise the sub-multiple algorithm tends
     * towards low F0 which leads to better sounding background noise.
     *
     * @param Sn a float array of input speech samples
     * @param prev_Wo_enc
     */
    public void nlp(float[] Sn, float prev_Wo_enc) {
        float notch;    // current notch filter output
        int i, j;

        /*
         * Square, notch filter at DC, and LP filter vector
         */
        for (i = (Defines.M - Defines.N); i < Defines.M; i++) {
            sq[i] = (Sn[i] * Sn[i]);
        }

        for (i = (Defines.M - Defines.N); i < Defines.M; i++) {     // notch filter at DC
            notch = sq[i] - mem_x;
            notch += COEFF * mem_y;
            mem_x = sq[i];
            mem_y = notch;

            /*
             * With 0 input vectors to codec, kiss_fft() would take a long
             * time to execute when running in real time.  Problem was traced
             * to kiss_fft function call in this function.
             * 
             * Adding this small constant fixed problem.  Not exactly sure why.
             */
            sq[i] = notch + 1.0F;
        }

        for (i = (Defines.M - Defines.N); i < Defines.M; i++) {     // FIR filter vector
            for (j = 0; j < (Nlp.NLP_NTAP - 1); j++) {
                mem_fir[j] = mem_fir[j + 1];
            }
            mem_fir[NLP_NTAP - 1] = sq[i];

            sq[i] = 0.0F;
            for (j = 0; j < Nlp.NLP_NTAP; j++) {
                sq[i] += mem_fir[j] * nlp_fir[j];
            }
        }

        // decimate and DFT
        for (i = 0; i < (Defines.M / Nlp.DEC); i++) {
            fw[i] = new Complex(sq[i * Nlp.DEC] * coswintab[i], 0.0F);
        }

        // fill the rest with 0.0
        for (i = (Defines.M / Nlp.DEC); i < Defines.FFT_SIZE; i++) {
            fw[i] = new Complex();
        }

        fftEncode.transform(fw);

        for (i = 0; i < Defines.FFT_SIZE; i++) {
            fw[i] = new Complex(fw[i].csqr(), fw[i].getImaginary());
        }

        /* find global peak */
        float gmax = 0.0F;
        int gmax_bin = Defines.FFT_SIZE * Nlp.DEC / Defines.P_MAX;      // 16

        /*
         * 512 * 5 / 160(P_MAX) = 16  and 512 * 5 / 20(P_MIN) = 128
         */
        for (i = Defines.FFT_SIZE * Nlp.DEC / Defines.P_MAX; i <= Defines.FFT_SIZE * Nlp.DEC / Defines.P_MIN; i++) {
            if (fw[i].getReal() > gmax) {
                gmax = fw[i].getReal();
                gmax_bin = i;
            }
        }

        /* Shift samples in buffer to make room for new samples */
        for (i = 0; i < (Defines.M - Defines.N); i++) {
            sq[i] = sq[i + Defines.N];
        }

        this.pitch = (float) (Defines.FS / post_process_sub_multiples(gmax, gmax_bin, prev_Wo_enc));
    }

    /**
     * Given the global maximma of Fw[] we search integer sub-multiples for local
     * maxima. If local maxima exist and they are above an experimentally
     * derived threshold (OK a magic number I pulled out of the air) we choose
     * the submultiple as the F0 estimate.
     *
     * <p>
     * The rational for this is that the lowest frequency peak of fw[] should be
     * F0, as fw[] can be considered the autocorrelation function of Sw[] (the
     * speech spectrum). However sometimes due to phase effects the lowest
     * frequency maxima may not be the global maxima.
     *
     * <p>
     * This works OK in practice and favors low F0 values in the presence of
     * background noise which means the sinusoidal codec does an OK job of
     * synthesizing the background noise. High F0 in background noise tends to
     * sound more periodic introducing annoying artifacts.
     *
     * @param gmax
     * @param gmax_bin
     * @param prev_Wo_enc
     * @return
     */
    private int post_process_sub_multiples(float gmax, int gmax_bin, float prev_Wo) {
        int cmax_bin, mult, b;
        int bmin, bmax, lmax_bin;
        float thresh, lmax;
        int prev_fo_bin = (int) (prev_Wo * Nlp.PREV_BIN);

        /*
         * post process estimate by searching submultiples
         */
        mult = 2;
        cmax_bin = gmax_bin;

        while ((gmax_bin / mult) >= Nlp.MIN_BIN) {
            b = (gmax_bin / mult);			// determine search interval
            bmin = (int) (0.8F * (float) b);
            bmax = (int) (1.2F * (float) b);

            if (bmin < Nlp.MIN_BIN) {
                bmin = Nlp.MIN_BIN;
            }

            /*
             * lower threshold to favour previous frames pitch estimate,
             * this is a form of pitch tracking
             */
            if ((prev_fo_bin > bmin) && (prev_fo_bin < bmax)) {
                thresh = Nlp.CNLP * 0.5F * gmax;
            } else {
                thresh = Nlp.CNLP * gmax;
            }

            lmax = 0.0F;
            lmax_bin = bmin;

            for (b = bmin; b <= bmax; b++) /* look for maximum in interval */ {
                if (fw[b].getReal() > lmax) {
                    lmax = fw[b].getReal();
                    lmax_bin = b;
                }
            }

            if (lmax > thresh) {
                if ((lmax > fw[lmax_bin - 1].getReal()) && (lmax > fw[lmax_bin + 1].getReal())) {
                    cmax_bin = lmax_bin;
                }
            }

            mult++;
        }

        return cmax_bin * Defines.FS / (Defines.FFT_SIZE * Nlp.DEC);    // best FO
    }
}
