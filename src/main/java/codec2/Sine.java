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
 * Sinusoidal analysis and synthesis functions
 *
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Sine {

    private static final float SIXTY = (float) (60.0 * Defines.TWO_PI / Defines.FS);
    private static final float V_THRESH = 6.0F;
    private static final float R = (float) (Defines.TWO_PI / Defines.FFT_SIZE);
    //
    private final Complex[] encodeSw;
    private final Complex[] decodeSw;
    //
    private final FFT fftEncode;
    private final FFT fftDecode;
    //
    private final Complex[] Sw;
    private final Complex[] Ew;

    public Sine(FFT ffte, FFT fftd) {
        this.fftEncode = ffte;
        this.fftDecode = fftd;
        this.encodeSw = new Complex[Defines.FFT_SIZE];
        this.decodeSw = new Complex[Defines.FFT_SIZE];
        this.Sw = new Complex[Defines.FFT_SIZE];
        this.Ew = new Complex[Defines.FFT_SIZE];
    }

    /**
     * Return the result of the decode FFT
     *
     * @return complex array of fft output
     */
    public Complex[] getDecodeDomain() {
        return this.decodeSw;
    }

    /**
     * Return the result of the encode FFT
     *
     * @return complex array of fft output
     */
    public Complex[] getEncodeDomain() {
        return this.encodeSw;
    }

    /**
     * Finds the DFT of the current speech input speech frame.
     *
     * @param Sn
     * @param w
     */
    public void dft_speech(float[] Sn, float[] w) {
        int i;

        /*
         * Initialize the working buffer to complex zero
         */
        for (i = 0; i < Defines.FFT_SIZE; i++) {
            this.encodeSw[i] = new Complex();
        }

        /*
         * Center analysis window on time axis, we need to arrange input
         * to FFT this way to make FFT phases correct
         */
        // move 2nd half to start of FFT input vector
        for (i = 0; i < Defines.NW / 2; i++) {
            this.encodeSw[i] = new Complex(Sn[i + Defines.M / 2] * w[i + Defines.M / 2], 0.0F);
        }

        // move 1st half to end of FFT input vector
        for (i = 0; i < Defines.NW / 2; i++) {
            this.encodeSw[Defines.FFT_SIZE - Defines.NW / 2 + i]
                    = new Complex(Sn[i + Defines.M / 2 - Defines.NW / 2] * w[i + Defines.M / 2 - Defines.NW / 2], 0.0F);
        }

        this.fftEncode.transform(encodeSw);       // now in frequency domain
    }

    /**
     * Refines the current pitch estimate using the harmonic sum pitch
     * estimation technique.
     *
     * @param model
     * @param Sw
     */
    public void two_stage_pitch_refinement(Model model, Complex[] Sw) {
        float pmin, pmax;	// pitch refinment minimum, maximum

        // Coarse refinement
        pmax = (float) Defines.TWO_PI / model.getWo() + 5.0F;
        pmin = (float) Defines.TWO_PI / model.getWo() - 5.0F;
        hs_pitch_refinement(model, Sw, pmin, pmax, 1.0F);    // step 1.0 course

        // Fine refinement
        pmax = (float) Defines.TWO_PI / model.getWo() + 1.0F;
        pmin = (float) Defines.TWO_PI / model.getWo() - 1.0F;
        hs_pitch_refinement(model, Sw, pmin, pmax, 0.25F);   // step 0.25 fine

        // Limit range
        if (model.getWo() < Defines.WO_MIN) {
            model.setWo(Defines.WO_MIN);
        } else if (model.getWo() > Defines.WO_MAX) {    // changed to if-else
            model.setWo(Defines.WO_MAX);
        }

        model.setL((int) Math.floor(Math.PI / (double) model.getWo()));     // no rounding ??
    }

    /**
     * Harmonic sum pitch refinement function.
     *
     * <p>
     * pmin pitch search range minimum pmax pitch search range maximum step
     * pitch search step size model current pitch estimate in model.Wo
     *
     * <p>
     * model refined pitch estimate in model.Wo
     *
     * @param model
     * @param pmin
     * @param pmax
     * @param pstep
     */
    private void hs_pitch_refinement(Model model, Complex[] Sw, float pmin, float pmax, float pstep) {
        int m;		// loop variable
        int b;		// bin for current harmonic center
        float E;	// energy for current pitch
        float Wo;	// current "test" fundamental freq.
        float Wom;	// Wo that maximizes E
        float Em;	// mamimum energy
        float p;	// current pitch

        // Initialization
        Wom = model.getWo();
        model.setL((int) ((float) Math.PI / Wom));	// use initial pitch est. for L
        Em = 0.0F;

        // Determine harmonic sum for a range of Wo values
        for (p = pmin; p <= pmax; p += pstep) {
            E = 0.0F;
            Wo = Defines.TWO_PI / p;

            // Sum harmonic magnitudes
            for (m = 1; m <= model.getL(); m++) {
                b = (int) ((float) m * Wo / Sine.R + 0.5);   // same a Wo * 1/R
                E += Sw[b].csqr();
            }

            // Compare to see if this is a maximum
            if (E > Em) {
                Em = E;
                Wom = Wo;
            }
        }

        model.setWo(Wom);
    }

    /**
     * Estimates the complex amplitudes of the harmonics.
     *
     * @param model
     * @param Sw
     * @param est_phase
     */
    public void estimate_amplitudes(Model model, Complex[] Sw, boolean est_phase) {
        int i, m;		// loop variables
        int am, bm;		// bounds of current harmonic
        int b;                  // DFT bin of centre of current harmonic
        float den;		// denominator of amplitude expression

        float tmp = model.getWo() / Sine.R;

        for (m = 1; m <= model.getL(); m++) {
            am = (int) (((float) m - 0.5F) * tmp + 0.5F);   // lower
            bm = (int) (((float) m + 0.5F) * tmp + 0.5F);   // upper
            b = (int) ((float) m * tmp + 0.5F);             // center

            // Estimate ampltude of harmonic
            den = 0.0F;

            for (i = am; i < bm; i++) {
                den += Sw[i].csqr();
            }

            model.setA(m, (float) Math.sqrt(den));

            if (est_phase == true) {

                /*
                 * Estimate phase of harmonic, this is expensive in CPU for
                 * embedded devices so we make it an option
                 */
                model.setPhi(m, (float) Math.atan2(Sw[b].getImaginary(), Sw[b].getReal()));
            }
        }
    }

    /**
     * Returns the error of the MBE cost function for a given F0.
     *
     * <p>
     * Note: I think a lot of the operations below can be simplified as
     * W[].getImaginary() = 0 and has been normalized such that den always
     * equals 1.
     *
     * @param model
     * @param S
     * @param W
     */
    public void est_voicing_mbe(Model model, Complex[] S, Complex[] W) {
        Complex Am;
        float den;
        int i, m, al, bl, offset;

        double signal = 1E-4;

        for (i = 1; i <= (model.getL() / 4); i++) {
            signal += (model.getA(i) * model.getA(i));
        }

        for (i = 0; i < Defines.FFT_SIZE; i++) {
            this.Sw[i] = new Complex();
            this.Ew[i] = new Complex();
        }

        double tmp = model.getWo() * Defines.FFT_SIZE / (float) Defines.TWO_PI;
        double error = 1E-4;

        /* Just test across the harmonics in the first 1000 Hz (L/4) */
        for (i = 1; i <= (model.getL() / 4); i++) {
            Am = new Complex();
            den = 0.0F;

            al = (int) Math.ceil(((double) i - 0.5) * tmp);     // no rounding ??
            bl = (int) Math.ceil(((double) i + 0.5) * tmp);

            /* Estimate amplitude of harmonic assuming harmonic is totally voiced */
            offset = (int) ((double) Defines.FFT_SIZE / 2.0 - (double) i * tmp + 0.5);

            for (m = al; m < bl; m++) {
                Am = Am.add(S[m].times(W[offset + m].getReal()));
                den += (W[offset + m].getReal() * W[offset + m].getReal());
            }

            Am = Am.divide(den);

            /*
             * Determine error between estimated harmonic and original
             */
            for (m = al; m < bl; m++) {
                this.Sw[m] = Am.times(W[offset + m].getReal());
                this.Ew[m] = S[m].minus(this.Sw[m]);
                error += Ew[m].csqr();
            }
        }

        if ((float) (10.0 * Math.log10(signal / error)) > Sine.V_THRESH) {
            model.setVoiced(true);
        } else {
            model.setVoiced(false);
        }

        /*
         * post processing, helps clean up some voicing errors
         *
         * Determine the ratio of low freqency to high frequency energy,
         * voiced speech tends to be dominated by low frequency energy,
         * unvoiced by high frequency. This measure can be used to
         * determine if we have made any gross errors.
         */
        double elow = 1E-4;
        double ehigh = 1E-4;

        for (i = 1; i <= (model.getL() / 2); i++) {
            elow += (model.getA(i) * model.getA(i));
        }

        for (i = (model.getL() / 2); i <= model.getL(); i++) {
            ehigh += (model.getA(i) * model.getA(i));
        }

        double eratio = 10.0 * Math.log10(elow / ehigh);

        if (model.getVoiced() == false) {
            /*
             * Look for Type 1 errors, strongly V speech that has been
             * accidentally declared UV
             */

            if (eratio > 10.0) {
                model.setVoiced(true);
            }
        } else if (model.getVoiced() == true) {
            /*
             * Look for Type 2 errors, strongly UV speech that has been
             * accidentally declared V
             */

            if (eratio < -10.0) {
                model.setVoiced(false);
            }

            /*
             * A common source of Type 2 errors is the pitch estimator
             * gives a low (50Hz) estimate for UV speech, which gives a
             * good match with noise due to the close harmoonic spacing.
             * These errors are much more common than people with 50Hz
             * pitch, so we have just a small eratio threshold.
             */
            if ((eratio < -4.0) && (model.getWo() <= Sine.SIXTY)) {
                model.setVoiced(false);
            }
        }
    }

    /**
     * Synthesize a speech signal in the frequency domain from the sinusoidal
     * model parameters. Uses overlap-add with a trapezoidal window to smoothly
     * interpolate between frames.
     *
     * @param model
     * @param Sn_
     * @param Pn
     */
    public void synthesise(Model model, float[] Sn_, float[] Pn) {
        int i, j, b;

        /* Update window by shifting 80 samples left (10 ms) */
        
        System.arraycopy(Sn_, Defines.N, Sn_, 0, Defines.N -1);    // 0 - 78
        Sn_[Defines.N -1] = 0.0F;

        for (i = 0; i < Defines.FFT_SIZE; i++) {
            this.decodeSw[i] = new Complex();
        }

        /*
         * Nov 2010 - found that synthesis using time domain cos() functions
         * gives better results for synthesis frames greater than 10ms.
         * 
         * Inverse FFT synthesis using a 512 pt FFT works well for 10ms window.
         * I think (but am not sure) that the problem is related to the
         * quantisation of the harmonic frequencies to the FFT bin size,
         * e.g. there is a  8000/512 Hz step between FFT bins.  For some reason
         * this makes the speech from longer frames > 10ms sound poor.
         * The effect can also be seen when synthesising test signals like
         * single sine waves, some sort of amplitude modulation at the
         * frame rate.
         * 
         * Another possibility is using a larger FFT size (1024 or 2048).
         */

        /*
         * Now set up frequency domain synthesized speech
         */
        for (i = 1; i <= model.getL(); i++) {
            b = (int) ((float) i * model.getWo() * (float) Defines.FFT_SIZE / Defines.TWO_PI + 0.5F);

            if (b > ((Defines.FFT_SIZE / 2) - 1)) {
                b = (Defines.FFT_SIZE / 2) - 1;
            }

            this.decodeSw[b] = new Complex(model.getA(i) * (float) Math.cos(model.getPhi(i)),
                    model.getA(i) * (float) Math.sin(model.getPhi(i)));
            this.decodeSw[Defines.FFT_SIZE - b] = this.decodeSw[b].conjugate();
        }

        // Perform inverse DFT
        this.fftDecode.itransform(this.decodeSw);

        // Overlap add to previous samples
        for (i = 0; i < (Defines.N - 1); i++) {
            Sn_[i] += this.decodeSw[Defines.FFT_SIZE - Defines.N + 1 + i].getReal() * Pn[i];
        }

        // put the new data on the end of the window
        for (i = Defines.N - 1, j = 0; i < (2 * Defines.N); i++, j++) {
            Sn_[i] = (this.decodeSw[j].getReal() * Pn[i]);
        }
    }
}
