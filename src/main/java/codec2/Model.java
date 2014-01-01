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

import java.util.Random;

/**
 * Speech object that contains the model of the sample frame
 *
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Model {

    private final static float R = Defines.TWO_PI / Defines.FFT_SIZE;
    private final static int MAX_HARMONIC = 80;        // maximum number of harmonics
    private final static float BG_THRESH = 40.0F;      // only consider low levels signals for bg_est
    private final static float BG_BETA = 0.1F;         // averaging filter constant
    private final static float BG_MARGIN = 6.0F;
    //
    private final FFT fft;
    private float Wo;                                  // fundamental frequency estimate in radians
    private int L;                                     // number of harmonics
    private final float[] A;                           // amplitiude of each harmonic
    private final float[] phi;                         // phase of each harmonic
    private boolean voiced;                            // if this frame is voiced

    /*
     * harmonics this far above BG noise are 
     * randomised.  Helped make bg noise less 
     * spikey (impulsive) for mmt1, but speech was
     * perhaps a little rougher.
     */
    private final Random rand;
    //
    private final Complex[] pw;
    private final Complex[] H;                         // LPC freq domain samples
    private final Complex[] Ex;                        // excitation samples
    private final Complex[] A_;                        // synthesised harmonic samples

    /*
     * Structure to hold model parameters for one 10 ms frame
     */
    public Model(FFT fftc) {
        this.fft = fftc;
        this.pw = new Complex[Defines.FFT_SIZE];
        this.H = new Complex[Model.MAX_HARMONIC + 1];
        this.Ex = new Complex[Model.MAX_HARMONIC + 1];
        this.phi = new float[Model.MAX_HARMONIC + 1];
        this.A = new float[Model.MAX_HARMONIC + 1];
        this.A_ = new Complex[Model.MAX_HARMONIC + 1];
        this.rand = new Random(System.currentTimeMillis());     // seed the noise generator
        this.Wo = 0.0F;
        this.L = 0;
        this.voiced = false;
    }

    public void reInit() {
        this.Wo = 0.0F;
        this.L = 0;
        this.voiced = false;

        for (int i = 0; i < (Model.MAX_HARMONIC + 1); i++) {
            this.A[i] = this.phi[i] = 0.0F;
        }
    }

    public boolean getVoiced() {
        return this.voiced;
    }

    public void setVoiced(boolean val) {
        this.voiced = val;
    }

    public float getA(int index) {
        return this.A[index];
    }

    public void setA(int index, float val) {
        this.A[index] = val;
    }

    public float getPhi(int index) {
        return this.phi[index];
    }

    public void setPhi(int index, float val) {
        this.phi[index] = val;
    }

    public int getL() {
        return this.L;
    }

    public void setL(int val) {
        this.L = val;
    }

    public float getWo() {
        return this.Wo;
    }

    public void setWo(float val) {
        this.Wo = val;
    }

    /**
     * Apply first harmonic LPC correction at decoder. This helps improve low
     * pitch males after LPC modeling, like hts1a and morig.
     */
    public void apply_lpc_correction() {
        if (this.Wo < ((float) Math.PI * 150.0F / 4000.0F)) {
            this.A[1] *= 0.032F;
        }
    }

    /**
     * Postfilter to improve sound quality for speech with high levels of
     * background noise. Unlike mixed-excitation models requires no bits to be
     * transmitted to handle background noise.
     *
     * <p>
     * The post filter is designed to help with speech corrupted by background
     * noise. The zero phase model tends to make speech with background noise
     * sound "clicky". With high levels of background noise the low level
     * inter-formant parts of the spectrum will contain noise rather than speech
     * harmonics, so modelling them as voiced (i.e. a continuous, non-random
     * phase track) is inaccurate.
     *
     * <p>
     * Some codecs (like MBE) have a mixed voicing model that breaks the
     * spectrum into voiced and unvoiced regions. Several bits/frame (5-12) are
     * required to transmit the frequency selective voicing information. Mixed
     * excitation also requires accurate voicing estimation (parameter
     * estimators always break occasionally under exceptional conditions).
     *
     * <p>
     * In our case we use a post filter approach which requires no additional
     * bits to be transmitted. The decoder measures the average level of the
     * background noise during unvoiced frames. If a harmonic is less than this
     * level it is made unvoiced by randomizing it's phases.
     *
     * <p>
     * This idea is rather experimental. Some potential problems that may
     * happen:
     *
     * <p>
     * 1/ If someone says "aaaaaaaahhhhhhhhh" will background estimator track up
     * to speech level? This would be a bad thing.
     *
     * <p>
     * 2/ If background noise suddenly disappears from the source speech does
     * estimate drop quickly? What is noise suddenly re-appears?
     *
     * <p>
     * 3/ Background noise with a non-flat spectrum. Current algorithm just
     * considers spectrum as a whole, but this could be broken up into bands,
     * each with their own estimator.
     *
     * <p>
     * 4/ Males and females with the same level of background noise. Check
     * performance the same. Changing Wo affects width of each band, may affect
     * bg energy estimates.
     *
     * <p>
     * 5/ Not sure what happens during long periods of voiced speech e.g.
     * "sshhhhhhh"
     * 
     * @param bg_est
     */
    public void postfilter(float[] bg_est) {
        int m;

        // determine average energy across spectrum
        float e = 1E-12F;

        for (m = 1; m <= this.L; m++) {
            e += (this.A[m] * this.A[m]);
        }

        e = (float) (10.0 * Math.log10(e / this.L));

        /*
         * If beneath threshold, update bg estimate.  The idea
         * of the threshold is to prevent updating during high level
         * speech.
         */
        if ((e < Model.BG_THRESH) && !this.voiced) {
            bg_est[0] = bg_est[0] * (1.0F - Model.BG_BETA) + (e * Model.BG_BETA);
        }

        /*
         * now mess with phases during voiced frames to make any harmonics
         * less then our background estimate unvoiced.
         */

        if (this.voiced == true) {
            for (m = 1; m <= this.L; m++) {
                if (this.A[m] < (float) (Math.pow(10.0, (bg_est[0] + Model.BG_MARGIN) / 20.0))) {
                    this.phi[m] =  (float) Defines.TWO_PI * this.rand.nextFloat(); // random value between 0.0 and TWO_PI
                }
            }
        }
    }

    /**
     * Samples the complex LPC synthesis filter spectrum at the harmonic
     * frequencies.
     *
     * @param aks
     */
    private void aks_to_H(float[] aks) {
        int i;

        // Determine DFT of A(exp(jw))
        for (i = 0; i <= Defines.LPC_ORD; i++) {
            this.pw[i] = new Complex(aks[i], 0.0F);
        }

        for (i = Defines.LPC_ORD + 1; i < Defines.FFT_SIZE; i++) {
            this.pw[i] = new Complex();
        }

        fft.transform(this.pw);

        // Sample magnitude and phase at harmonics
        this.H[0] = new Complex();   // not used, but needs allocating

        float tmp = this.Wo / Model.R;

        for (int m = 1; m <= this.L; m++) {
            int am = (int) (((float) m - 0.5F) * tmp + 0.5F); // these used to have a floor()
            int bm = (int) (((float) m + 0.5F) * tmp + 0.5F);
            int b = (int) ((float) m * tmp + 0.5F);

            float Em = 0.0F;

            for (i = am; i < bm; i++) {
                Em += 1.0F / pw[i].csqr();     // was G / pw, but G was always one
            }

            float Am = (float) Math.sqrt(Math.abs(Em / (bm - am)));
            float phi_ = (float) -Math.atan2(this.pw[b].getImaginary(), this.pw[b].getReal());

            this.H[m] = new Complex((Am * (float) Math.cos(phi_)), (Am * (float) Math.sin(phi_)));
        }
    }

    /**
     * Synthesizes phases based on SNR and a rule based approach. No phase
     * parameters are required apart from the SNR (which can be reduced to a 1
     * bit V/UV decision per frame).
     *
     * <p>
     * The phase of each harmonic is modeled as the phase of a LPC synthesis
     * filter excited by an impulse. Unlike the first order model the position
     * of the impulse is not transmitted, so we create an excitation pulse train
     * using a rule based approach. * Consider a pulse train with a pulse
     * starting time n=0, with pulses repeated at a rate of Wo, the fundamental
     * frequency. A pulse train in the time domain is equivalent to harmonics in
     * the frequency domain. We can make an excitation pulse train using a sum
     * of sinsusoids:
     *
     * <p>
     * for (m=1; m <= L; m++)<br> ex[n] = cos(m * Wo * n)
     *
     * <p>
     * Note: the Octave script ../octave/phase.m is an example of this if you
     * would like to try making a pulse train.
     *
     * <p>
     * The phase of each excitation harmonic is:
     *
     * <p>
     * arg(E[m]) = mWo
     *
     * <p>
     * where E[m] are the complex excitation (freq domain) samples, arg(x), just
     * returns the phase of a complex sample x.
     *
     * <p>
     * As we don't transmit the pulse position for this model, we need to
     * synthesize it. Now the excitation pulses occur at a rate of Wo. This
     * means the phase of the first harmonic advances by N samples over a
     * synthesis frame of N samples. For example if Wo is pi/20 (200 Hz), then
     * over a 10ms frame (N=80 samples), the phase of the first harmonic would
     * advance (pi/20)*80 = 4*pi or two complete cycles.
     *
     * <p>
     * We generate the excitation phase of the fundamental (first harmonic):
     *
     * <p>
     * arg[E[1]] = Wo*N;
     *
     * <p>
     * We then relate the phase of the m-th excitation harmonic to the phase of
     * the fundamental as:
     *
     * <p>
     * arg(E[m]) = m*arg(E[1])
     *
     * <p>
     * This E[m] then gets passed through the LPC synthesis filter to determine
     * the final harmonic phase.
     *
     * <p>
     * Comparing to speech synthesized using original phases:
     *
     * <p>
     * - Through headphones speech synthesized with this model is not as good.
     * Through a loudspeaker it is very close to original phases.
     *
     * <p>
     * - If there are voicing errors, the speech can sound clicky or staticy. If
     * V speech is mistakenly declared UV, this model tends to synthesize
     * impulses or clicks, as there is usually very little shift or dispersion
     * through the LPC filter.
     *
     * <p>
     * - When combined with LPC amplitude modeling there is an additional drop
     * in quality. I am not sure why, theory is interformant energy is raised
     * making any phase errors more obvious.
     *
     * <p>
     * NOTES:
     *
     * <p>
     * 1/ This synthesis model is effectively the same as a simple LPC-10
     * vocoders, and yet sounds much better. Why? Conventional wisdom (AMBE,
     * MELP) says mixed voicing is required for high quality speech.
     *
     * <p>
     * 2/ I am pretty sure the Lincoln Lab sinusoidal coding guys (like xMBE
     * also from MIT) first described this zero phase model, I need to look up
     * the paper.
     *
     * <p>
     * 3/ Note that this approach could cause some discontinuities in the phase
     * at the edge of synthesis frames, as no attempt is made to make sure that
     * the phase tracks are continuous (the excitation phases are continuous,
     * but not the final phases after filtering by the LPC spectra). Technically
     * this is a bad thing. However this may actually be a good thing,
     * disturbing the phase tracks a bit. More research needed, e.g. test a
     * synthesis model that adds a small delta-W to make phase tracks line up
     * for voiced harmonics.
     *
     * @param aks
     * @param ex_phase
     */
    public void phase_synth_zero_order(float[] aks, float[] ex_phase) {
        aks_to_H(aks);

        /*
         * Update excitation fundamental phase track, this sets the position
         * of each pitch pulse during voiced speech.  After much experiment
         * I found that using just this frame's Wo improved quality for UV
         * sounds compared to interpolating two frames Wo like this:
         * 
         * ex_phase[0] += (*prev_Wo+model->Wo)*N/2;
         */
        ex_phase[0] += (this.Wo * Defines.N);
        ex_phase[0] -= (float) (Defines.TWO_PI * Math.floor(ex_phase[0] / Defines.TWO_PI + 0.5));

        this.A_[0] = new Complex(); // allocate, but not used
        this.Ex[0] = new Complex(); // ditto

        for (int m = 1; m <= this.L; m++) {

            /*
             * generate excitation
             */
            if (this.voiced == true) {
                this.Ex[m] = new Complex((float) Math.cos(ex_phase[0] * (double) m), (float) Math.sin(ex_phase[0] * (double) m));
            } else {

                /*
                 * When a few samples were tested I found that LPC filter
                 * phase is not needed in the unvoiced case, but no harm in
                 * keeping it.
                 */
                float lphi = (float) Defines.TWO_PI * this.rand.nextFloat();  // random value between 0.0 and TWO_PI

                this.Ex[m] = new Complex((float) Math.cos(lphi), (float) Math.sin(lphi));
            }

            /*
             * filter using LPC filter
             */
            this.A_[m] = this.H[m].times(this.Ex[m]);

            /*
             * modify sinusoidal phase
             */
            this.phi[m] = (float) Math.atan2(A_[m].getImaginary(), A_[m].getReal() + 1E-12);
        }
    }
}
