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
 * Quantization functions for the sinusoidal coder
 *
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Quantise {

    private static final float E_MIN_DB = -10.0F;
    private static final float E_MAX_DB = 40.0F;
    private static final float E_MAX_MINUS_E_MIN = (E_MAX_DB - E_MIN_DB);
    private static final float LSP_DELTA = 0.01F;                               // grid spacing for LSP root searches
    private static final float R = ((float) Defines.TWO_PI / Defines.FFT_SIZE); // number of rads/bin
    private static final float ESTEP = (E_MAX_MINUS_E_MIN / Defines.E_LEVELS);
    private static final float OL = (float) (50.0 * (Math.PI / 4000.0));
    private static final float OH = (float) (100.0 * (Math.PI / 4000.0));
    //
    private final FFT fftDecode;
    private final Complex[] pw;
    private final Complex[] xa;                                                // input to FFTs
    private final float[] Rw;                                                  // R = WA
    private final float[] ge_coeff = new float[2];
    private final Lsp lspclass;

    public Quantise(FFT fft, Lsp lsp) {
        fftDecode = fft;
        this.pw = new Complex[Defines.FFT_SIZE];
        this.xa = new Complex[Defines.FFT_SIZE];
        this.Rw = new float[Defines.FFT_SIZE];
        this.ge_coeff[0] = 0.8F;
        this.ge_coeff[1] = 0.9F;

        this.lspclass = lsp;
    }

    /**
     *
     * @param i
     * @return
     */
    public int lsp_bits(int i) {
        return Codebookcb.lsp_cb[i].log2m;
    }

    /**
     *
     * @param i
     * @return
     */
    public int lspd_bits(int i) {
        return Codebookd.lsp_cbd[i].log2m;
    }

    /**
     *
     * @param i
     * @return
     */
    public int lsp_pred_vq_bits(int i) {
        return Codebookjvm.lsp_cbjvm[i].log2m;
    }

    /**
     * Thirty-six bit scalar LSP quantizer. From a vector of un-quantized LSPs
     * finds the quantized LSP indexes.
     *
     * @param indexes
     * @param lsp
     */
    public void encode_lsps_scalar(int[] indexes, float[] lsp) {
        float[] lsp_hz = new float[Defines.LPC_ORD];
        float[] cb;
        int i, m;

        /*
         * convert from radians to Hz so we can use human readable
         * frequencies
         */
        for (i = 0; i < Defines.LPC_ORD; i++) {
            lsp_hz[i] =  (4000.0F / (float) Math.PI) * lsp[i];
        }

        /* scalar quantisers */
        for (i = 0; i < Defines.LPC_ORD; i++) {
            m = Codebookcb.lsp_cb[i].m;
            cb = Codebookcb.lsp_cb[i].cb;
            indexes[i] = quantise(cb, lsp_hz[i], m);
        }
    }

    /**
     * From a vector of quantized LSP indexes, returns the quantized LSPs.
     *
     * @param lsp
     * @param indexes
     */
    public void decode_lsps_scalar(float[] lsp, int[] indexes) {
        float[] lsp_hz = new float[Defines.LPC_ORD];
        float[] cb;
        int i, k;

        for (i = 0; i < Defines.LPC_ORD; i++) {
            k = Codebookcb.lsp_cb[i].k;
            cb = Codebookcb.lsp_cb[i].cb;
            lsp_hz[i] = cb[indexes[i] * k];
        }

        // convert back to radians
        for (i = 0; i < Defines.LPC_ORD; i++) {
            lsp[i] = ((float)Math.PI / 4000.0F) * lsp_hz[i];
        }

        sort_lsp_order(lsp);
        bw_expand_lsps(lsp);
    }

    /**
     * Scalar/VQ LSP difference quantizer.
     *
     * @param indexes
     * @param lsp
     */
    public void encode_lspds_scalar(int[] indexes, float[] lsp) {
        float[] lsp_hz = new float[Defines.LPC_ORD];
        float[] lsp__hz = new float[Defines.LPC_ORD];
        float[] dlsp = new float[Defines.LPC_ORD];
        float[] dlsp_ = new float[Defines.LPC_ORD];
        float[] cb;
        int i, m;

        /*
         * convert from radians to Hz so we can use human readable
         * frequencies
         */
        for (i = 0; i < Defines.LPC_ORD; i++) {
            lsp_hz[i] =  (4000.0F / (float) Math.PI) * lsp[i];
        }

        for (i = 0; i < Defines.LPC_ORD; i++) {

            /*
             * find difference from previous qunatised lsp
             */
            if (i != 0) {
                dlsp[i] = lsp_hz[i] - lsp__hz[i - 1];
            } else {
                dlsp[0] = lsp_hz[0];
            }

            m = Codebookd.lsp_cbd[i].m;
            cb = Codebookd.lsp_cbd[i].cb;
            indexes[i] = quantise(cb, dlsp[i], m);

            dlsp_[i] = cb[indexes[i]];

            if (i != 0) {
                lsp__hz[i] = lsp__hz[i - 1] + dlsp_[i];
            } else {
                lsp__hz[0] = dlsp_[0];
            }
        }
    }

    /**
     *
     * @param lsp
     * @param indexes
     */
    public void decode_lspds_scalar(float[] lsp, int[] indexes) {
        float[] cb;
        float[] lsp_hz = new float[Defines.LPC_ORD];
        float[] dlsp = new float[Defines.LPC_ORD];
        int i;

        for (i = 0; i < Defines.LPC_ORD; i++) {
            cb = Codebookd.lsp_cbd[i].cb;
            dlsp[i] = cb[indexes[i]];

            if (i != 0) {
                lsp_hz[i] = lsp_hz[i - 1] + dlsp[i];
            } else {
                lsp_hz[0] = dlsp[0];
            }

            lsp[i] = ((float) Math.PI / 4000.0F) * lsp_hz[i];
        }
    }

    /**
     *
     * @param codebook
     * @param nb_entries
     * @param x
     * @return
     */
    private int find_nearest(float[] codebook, int nb_entries, float[] x) {
        int i, j, v;
        float val, dist, min_dist = 1e15F;
        int nearest = 0;

        for (i = 0; i < nb_entries; i++) {
            dist = 0.0F;

            for (j = 0; j < Defines.LPC_ORD; j++) {
                v = i * Defines.LPC_ORD + j;
                val = x[j] - codebook[v];
                dist += (val * val);
            }

            if (dist < min_dist) {
                min_dist = dist;
                nearest = i;
            }
        }

        return nearest;
    }

    /**
     *
     * @param codebook
     * @param nb_entries
     * @param err
     * @param w
     * @param ndim
     * @return
     */
    private int find_nearest_weighted(float[] codebook, int nb_entries, float[] err, float[] w, int ndim) {
        float val, dist, min_dist = 1e15F;
        int nearest = 0;

        for (int i = 0; i < nb_entries; i++) {
            dist = 0.0F;

            for (int j = 0; j < ndim; j++) {
                val = err[j] - codebook[i * ndim + j];
                dist += (w[j] * (val * val));
            }

            if (dist < min_dist) {
                min_dist = dist;
                nearest = i;
            }
        }

        return nearest;
    }

    /**
     * Applies a post filter to the LPC synthesis filter power spectrum Pw,
     * which suppresses the inter-formant energy.
     *
     * <p>
     * The algorithm is from p267 (Section 8.6) of "Digital Speech", edited by
     * A.M. Kondoz, 1994 published by Wiley and Sons. Chapter 8 of this text is
     * on the MBE vocoder, and this is a freq domain adaptation of post
     * filtering commonly used in CELP.
     *
     * <p>
     * I used the Octave simulation lpcpf.m to get an understanding of the
     * algorithm.
     *
     * <p>
     * Requires two more FFTs which is significantly more MIPs. However it
     * should be possible to implement this more efficiently in the time domain.
     * Just not sure how to handle relative time delays between the synthesis
     * stage and updating these coeffs. A smaller FFT size might also be
     * acceptable to save CPU.
     *
     * <p>
     * TODO: [ ] sync var names between Octave and C version [ ] doc gain
     * normalization [ ] I think the first FFT is not rqd as we do the same
     * thing in aks_to_M2().
     *
     * @param ak
     */
    private void lpc_post_filter(float[] ak, boolean bass_boost, float beta, float gamma) {
        int i;
        float e_before, e_after, gain;
        float Pfw;
        float[] Wmag = new float[Defines.FFT_SIZE / 2];
        float[] Amag = new float[Defines.FFT_SIZE / 2];
        float max_Rw, min_Rw;

        /*
         * Determine LPC inverse filter spectrum 1/A(exp(jw))
         *
         * We actually want the synthesis filter A(exp(jw)) but the
         * inverse (analysis) filter is easier to find as it's FIR, we
         * just use the inverse of 1/A to get the synthesis filter
         * A(exp(jw))
         */
        for (i = 0; i <= Defines.LPC_ORD; i++) {
            this.xa[i] = new Complex(ak[i], 0.0F);
        }

        for (i = Defines.LPC_ORD + 1; i < Defines.FFT_SIZE; i++) {
            this.xa[i] = new Complex();
        }

        fftDecode.transform(xa);

        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            Amag[i] = 1.0F / xa[i].csqr();
        }

        // Determine weighting filter spectrum W(exp(jw))
        this.xa[0] = new Complex(ak[0], 0.0F);
        float coef = gamma;

        for (i = 1; i <= Defines.LPC_ORD; i++) {
            this.xa[i] = new Complex(ak[i] * coef, 0.0F);
            coef *= gamma;
        }

        for (i = Defines.LPC_ORD + 1; i < Defines.FFT_SIZE; i++) {
            this.xa[i] = new Complex();
        }

        fftDecode.transform(this.xa);

        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            Wmag[i] = xa[i].csqr();
        }

        // Determined combined filter R = WA
        max_Rw = 0.0F;
        min_Rw = 1E32F;

        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            this.Rw[i] = (float) Math.sqrt(Wmag[i] * Amag[i]);

            if (this.Rw[i] > max_Rw) {
                max_Rw = Rw[i];
            }

            if (this.Rw[i] < min_Rw) {
                min_Rw = this.Rw[i];
            }
        }

        // create post filter mag spectrum and apply
        // measure energy before post filtering
        e_before = 1E-4F;
        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            e_before += this.pw[i].getReal();
        }

        // apply post filter and measure energy
        e_after = 1E-4F;

        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            Pfw =  (float) Math.pow(this.Rw[i], beta);
            this.pw[i] = this.pw[i].times(Pfw * Pfw);
            e_after += this.pw[i].getReal();
        }

        gain = e_before / e_after;

        // apply gain factor to normalise energy
        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            this.pw[i] = new Complex(this.pw[i].getReal() * gain, this.pw[i].getImaginary());
        }

        if (bass_boost == true) {
            // add 3dB to first 1 kHz to account for LP effect of PF

            for (i = 0; i < Defines.FFT_SIZE / 8; i++) {
                this.pw[i] = new Complex(this.pw[i].getReal() * (1.4F * 1.4F), this.pw[i].getImaginary());
            }
        }
    }

    /**
     * Transforms the linear prediction coefficients to spectral amplitude
     * samples. This function determines A(m) from the average energy per band
     * using an FFT.
     *
     * @param ak
     * @param model
     * @param E
     * @param lpc_pf
     * @param bass_boost
     * @param beta
     * @param gamma
     * @param snr a float representing the signal to noise ratio
     */
    public void aks_to_M2(float[] ak, Model model, float E, boolean lpc_pf, boolean bass_boost,
            float beta, float gamma, float[] snr) {
        int am, bm;		// limits of current band
        float Em;		// energy in band
        int i, m;

        // Determine DFT of A(exp(jw))
        for (i = 0; i <= Defines.LPC_ORD; i++) {
            this.pw[i] = new Complex(ak[i], 0.0F);
        }

        for (i = Defines.LPC_ORD + 1; i < Defines.FFT_SIZE; i++) {
            this.pw[i] = new Complex();
        }

        fftDecode.transform(pw);

        // Determine power spectrum P(w) = E/(A(exp(jw))^2
        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            this.pw[i] = new Complex(E / pw[i].csqr(), pw[i].getImaginary());
        }

        if (lpc_pf == true) {
            lpc_post_filter(ak, bass_boost, beta, gamma);   // result in static pw[]
        }

        // Determine magnitudes from P(w)
        /*
         * when used just by decoder {A} might be all zeroes so init signal
         * and noise to prevent log(0) errors
         */
        double signal = 1E-30F;
        double noise = 1E-32F;
        float tmp = model.getWo() / R;

        for (m = 1; m <= model.getL(); m++) {
            am = (int) (((float) m - 0.5F) * tmp + 0.5F);
            bm = (int) (((float) m + 0.5F) * tmp + 0.5F);
            Em = 0.0F;

            for (i = am; i < bm; i++) {
                Em += this.pw[i].getReal();
            }

            float Am = (float) Math.sqrt(Em);

            signal += (model.getA(m) * model.getA(m));
            noise += (model.getA(m) - Am) * (model.getA(m) - Am);

            /*
             * This code significantly improves perf of LPC model, in
             * particular when combined with phase0.  The LPC spectrum tends
             * to track just under the peaks of the spectral envelope, and
             * just above nulls.  This algorithm does the reverse to
             * compensate - raising the amplitudes of spectral peaks, while
             * attenuating the null.  This enhances the formants, and
             * supresses the energy between formants.
             */
            model.setA(m, Am);
        }

        snr[0] = (float) (10.0 * Math.log10(signal / noise));
    }

    /**
     * Encodes Wo using a WO_LEVELS quantizer.
     *
     * @param Wo
     * @return
     */
    public int encode_Wo(float Wo) {
        double norm = (Wo - Defines.WO_MIN) / (Defines.WO_MAX - Defines.WO_MIN);
        int index = (int) Math.floor((double) Defines.WO_LEVELS * norm + 0.5);

        if (index < 0) {
            index = 0;
        } else if (index > (Defines.WO_LEVELS - 1)) {
            index = Defines.WO_LEVELS - 1;
        }

        return index;
    }

    /*
     * Decodes Wo using a WO_LEVELS quantiser.
     */
    public float decode_Wo(int index) {
        return (Defines.WO_MIN + (Defines.STEP * (float) index));
    }

    /**
     * Analyze a windowed frame of time domain speech to determine LPCs which
     * are the converted to LSPs for quantization and transmission over the
     * channel.
     *
     * @param lsp
     * @param Sn
     * @param w
     * @return
     */
    public float speech_to_uq_lsps(float[] lsp, float[] Sn, float[] w) {
        float[] Wn = new float[Defines.M];
        float[] ak;
        float[] R1;
        float e;
        int i, roots;

        e = 0.0F;
        for (i = 0; i < Defines.M; i++) {
            Wn[i] = Sn[i] * w[i];
            e += (Wn[i] * Wn[i]);
        }

        /* trap 0 energy case as LPC analysis will fail */
        if (e == 0.0) {
            for (i = 0; i < Defines.LPC_ORD; i++) {
                lsp[i] = ((float) Math.PI / Defines.LPC_ORD) * (float) i;
            }

            return 0.0F;
        }

        R1 = autocorrelate(Wn);
        ak = levinson_durbin(R1);

        e = 0.0F;
        for (i = 0; i <= Defines.LPC_ORD; i++) {
            e += (ak[i] * R1[i]);
        }

        /*
         * 15 Hz BW expansion as I can't hear the difference and it may help
         * help occasional fails in the LSP root finding.  Important to do this
         * after energy calculation to avoid -ve energy values.
         */
        for (i = 0; i <= Defines.LPC_ORD; i++) {
            ak[i] *= (float) Math.pow(0.994, (double) i);
        }

        roots = lspclass.lpc_to_lsp(ak, lsp, 5, Quantise.LSP_DELTA);

        if (roots != Defines.LPC_ORD) {
            /* if root finding fails use some benign LSP values instead */
            for (i = 0; i < Defines.LPC_ORD; i++) {
                lsp[i] = ((float) Math.PI / (float) Defines.LPC_ORD) * (float) i;
            }
        }

        return e;
    }

    /**
     * Quantizes vector by choosing the nearest vector in codebook, and returns
     * the vector index.
     *
     * @param cb current VQ codebook
     * @param vec vector to quantize
     * @param m size of codebook
     * @return an int representing the best index so far
     */
    private int quantise(float[] cb, float vec, int m) {
        float e;	// current error
        int best_i;	// best index so far
        float best_e;	// best error so far
        float diff;

        best_i = 0;
        best_e = 1E32F;

        for (int j = 0; j < m; j++) {
            diff = cb[j] - vec;
            e = (diff * diff);       // pow(x,2) is slow

            if (e < best_e) {
                best_e = e;
                best_i = j;
            }
        }

        return best_i;
    }

    /**
     * Encodes LPC energy using an E_LEVELS quantizer.
     *
     * @param e
     * @return
     */
    public int encode_energy(float e) {
        double norm = ((10.0 * Math.log10((double) e)) - E_MIN_DB) / E_MAX_MINUS_E_MIN;
        int index = (int) Math.floor(Defines.E_LEVELS * norm + 0.5);

        if (index < 0) {
            index = 0;
        } else if (index > (Defines.E_LEVELS - 1)) {
            index = (Defines.E_LEVELS - 1);
        }

        return index;
    }

    /**
     * Decodes energy using a E_LEVELS quantizer.
     *
     * @param index
     * @return
     */
    public float decode_energy(int index) {
        return (float) Math.pow(10.0, (double) (E_MIN_DB + ESTEP * (float) index) / 10.0F);
    }

    /**
     *
     * @param xl
     * @param w
     */
    private void compute_weights(float[] xl, float[] w, float[] xq) {
        w[0] = 30.0F;
        w[1] = 1.0F;

        if (xl[1] < 0.0) {
            w[0] *= .6;
            w[1] *= .3;
        } else if (xl[1] < -10.0) {
            w[0] *= .3;
            w[1] *= .3;
        }

        /* Higher weight if pitch is stable */
        if (Math.abs(xl[0] - xq[0]) < .2) {
            w[0] *= 2.0;
            w[1] *= 1.5;
        } else if (Math.abs(xl[0] - xq[0]) > .5) {          // Lower if not stable
            w[0] *= .5;
        }

        /* Lower weight for low energy */
        if (xl[1] < xq[1] - 10.0) {
            w[1] *= .5;
        } else if (xl[1] < xq[1] - 20.0) {
            w[1] *= .5;
        }

        /* Square the weights because it's applied on the squared error */
        w[0] *= w[0];
        w[1] *= w[1];
    }

    /**
     * Joint Wo and LPC energy vector quantizer developed by Jean-Marc Valin.
     * Returns index, and updated states xq[].
     *
     * @param model
     * @param e
     * @param err
     * @param xq
     * @return
     */
    public int encode_WoE(Model model, float e, float[] err, float[] xq) {
        int i, n1;
        int nb_entries = Codebookge.lsp_cbge[0].m;      // 256
        int ndim = Codebookge.lsp_cbge[0].k;            // 2
        float[] codebook1 = Codebookge.lsp_cbge[0].cb;
        float[] xl = new float[ndim];
        float[] w = new float[ndim];

        if (e < 0.0) {
            e = 0.0F;  // occasional small negative energies due LPC round off I guess
        }

        xl[0] = (float) (Math.log10((model.getWo() / Math.PI) * 4000.0 / 50.0) / Math.log10(2.0));
        xl[1] = (float) (10.0 * Math.log10(1e-4 + e));

        compute_weights(xl, w, xq);

        for (i = 0; i < ndim; i++) {
            err[i] = xl[i] - this.ge_coeff[i] * xq[i];
        }

        n1 = find_nearest_weighted(codebook1, nb_entries, err, w, ndim);

        for (i = 0; i < ndim; i++) {
            xq[i] = this.ge_coeff[i] * xq[i] + codebook1[ndim * n1 + i];
            err[i] -= codebook1[ndim * n1 + i];
        }

        return n1;
    }

    /**
     * Joint Wo and LPC energy vector quantizer developed by Jean-Marc Valin.
     * Given index and states xq[], returns Wo & E, and updates states xq[].
     *
     * @param model
     * @param n1
     * @param xq
     * @return
     */
    public float decode_WoE(Model model, int n1, float[] xq) {
        int ndim = Codebookge.lsp_cbge[0].k;

        for (int i = 0; i < ndim; i++) {
            xq[i] = this.ge_coeff[i] * xq[i] + Codebookge.lsp_cbge[0].cb[ndim * n1 + i];
        }

        float Wo = (float) (Math.pow(2.0, xq[0]) * (Math.PI * 50.0) / 4000.0);

        /*
         * bit errors can make us go out of range leading to all sorts of
         * probs like segmentation faults
         */
        if (Wo > Defines.WO_MAX) {
            Wo = Defines.WO_MAX;
        } else if (Wo < Defines.WO_MIN) {
            Wo = Defines.WO_MIN;
        }

        model.setL((int) ((float) Math.PI / Wo)); // re-compute L
        model.setWo(Wo);

        return ((float) Math.pow(10.0, xq[1] / 10.0));
    }

    /**
     * Multi-stage VQ LSP quantizer developed by Jean-Marc Valin.
     *
     * @param indexes
     * @param lsp
     */
    public void encode_lsps_vq(int[] indexes, float[] lsp) {
        float[] err = new float[Defines.LPC_ORD];
        float[] err2 = new float[Defines.LPC_ORD / 2];
        float[] err3 = new float[Defines.LPC_ORD / 2];
        float[] w = new float[Defines.LPC_ORD];
        float[] w2 = new float[Defines.LPC_ORD / 2];
        float[] w3 = new float[Defines.LPC_ORD / 2];
        float[] codebook1 = Codebookjvm.lsp_cbjvm[0].cb;
        float[] codebook2 = Codebookjvm.lsp_cbjvm[1].cb;
        float[] codebook3 = Codebookjvm.lsp_cbjvm[2].cb;
        int i;

        w[0] = Math.min(lsp[0], lsp[1] - lsp[0]);

        for (i = 1; i < (Defines.LPC_ORD - 1); i++) {
            w[i] = Math.min(lsp[i] - lsp[i - 1], lsp[i + 1] - lsp[i]);
        }

        w[Defines.LPC_ORD - 1] = Math.min(lsp[Defines.LPC_ORD - 1] - lsp[Defines.LPC_ORD - 2], (float) Math.PI - lsp[Defines.LPC_ORD - 1]);

        for (i = 0; i < Defines.LPC_ORD; i++) {
            w[i] = 1.0F / (.01F + w[i]);
        }

        int n1 = find_nearest(codebook1, Codebookjvm.lsp_cbjvm[0].m, lsp);

        for (i = 0; i < Defines.LPC_ORD; i++) {
            err[i] = lsp[i] - codebook1[Defines.LPC_ORD * n1 + i];
        }

        for (i = 0; i < (Defines.LPC_ORD / 2); i++) {
            err2[i] = err[2 * i];
            err3[i] = err[2 * i + 1];
            w2[i] = w[2 * i];
            w3[i] = w[2 * i + 1];
        }

        int n2 = find_nearest_weighted(codebook2, Codebookjvm.lsp_cbjvm[1].m, err2, w2, Defines.LPC_ORD / 2);
        int n3 = find_nearest_weighted(codebook3, Codebookjvm.lsp_cbjvm[2].m, err3, w3, Defines.LPC_ORD / 2);

        indexes[0] = n1;
        indexes[1] = n2;
        indexes[2] = n3;
    }

    /**
     *
     * @param lsp
     * @param indexes
     */
    public void decode_lsps_vq(float[] lsp, int[] indexes) {
        float[] codebook1 = Codebookjvm.lsp_cbjvm[0].cb;
        float[] codebook2 = Codebookjvm.lsp_cbjvm[1].cb;
        float[] codebook3 = Codebookjvm.lsp_cbjvm[2].cb;
        int n1 = indexes[0];
        int n2 = indexes[1];
        int n3 = indexes[2];
        int i;

        for (i = 0; i < Defines.LPC_ORD; i++) {
            lsp[i] = codebook1[Defines.LPC_ORD * n1 + i];
        }

        for (i = 0; i < (Defines.LPC_ORD / 2); i++) {
            lsp[2 * i] += codebook2[Defines.LPC_ORD * n2 / 2 + i];
            lsp[2 * i + 1] += codebook3[Defines.LPC_ORD * n3 / 2 + i];
        }

        sort_lsp_order(lsp);
        bw_expand_lsps(lsp);
    }

    private void sort_lsp_order(float[] lsp) {
        boolean swapped;
        float tmp;

        do {
            swapped = false;

            for (int i = 1; i < Defines.LPC_ORD; i++) {
                if (lsp[i] < lsp[i - 1]) {
                    tmp = lsp[i - 1];               // swap i and i-1
                    lsp[i - 1] = lsp[i] - 0.1F;
                    lsp[i] = tmp + 0.1F;
                    swapped = true;
                }
            }
        } while (swapped);      // sooner or later they should all be sorted
    }

    private void bw_expand_lsps(float[] lsp) {
        int i;

        /*
         * Apply BW Expansion to this sorted vector of LSPs.
         */
        for (i = 1; i < Defines.LPC_ORD; i++) {

            /*
             * Prevents any two LSPs getting too close together after quantisation.
             * We know from experiment that LSP quantisation errors < 12.5Hz (25Hz
             * step size) are inaudible so we use that as the minimum LSP separation.
             */
            if (i < 4) {
                if ((lsp[i] - lsp[i - 1]) < Quantise.OL) {    // 50x
                    lsp[i] = lsp[i - 1] + Quantise.OL;
                }
            } else {

                /*
                 * As quantiser gaps increase, larger BW expansion was required
                 * to prevent twinkly noises.  This may need more experiment for
                 * different quanstisers.
                 */
                if (lsp[i] - lsp[i - 1] < Quantise.OH) {     // 100x
                    lsp[i] = lsp[i - 1] + Quantise.OH;
                }
            }
        }
    }

    /**
     * Finds the first P autocorrelation values of an array of windowed speech
     * samples.
     *
     * @param Sn a float array of windowed speech samples
     * @return a float array of P+1 autocorrelation coefficients
     */
    private float[] autocorrelate(float[] Sn) {
        float[] Rn = new float[Defines.LPC_ORD + 1];

        for (int j = 0; j < Defines.LPC_ORD + 1; j++) {
            Rn[j] = 0.0F;

            for (int i = 0; i < (Defines.M - j); i++) {
                Rn[j] += (Sn[i] * Sn[i + j]);
            }
        }

        return Rn;
    }

    /**
     * Given P+1 autocorrelation coefficients, finds P Linear Prediction Coeff.
     * (LPCs) where P is the order of the LPC all-pole model. The
     * Levinson-Durbin algorithm is used, and is described in:
     *
     * J. Makhoul "Linear prediction, a tutorial review" Proceedings of the IEEE
     * Vol-63, No. 4, April 1975
     *
     * @param r
     * @return
     */
    private float[] levinson_durbin(float[] r) {
        float[][] a = new float[Defines.LPC_ORD + 1][Defines.LPC_ORD + 1];
        float[] ak = new float[Defines.LPC_ORD + 1];
        float[] E = new float[Defines.LPC_ORD + 1];
        float[] k = new float[Defines.LPC_ORD + 1];
        float sum;
        int i, j;

        E[0] = r[0];                                            // Equation 38a, Makhoul

        for (i = 1; i <= Defines.LPC_ORD; i++) {
            sum = 0.0F;

            for (j = 1; j <= i - 1; j++) {
                sum += (a[i - 1][j] * r[i - j]);
            }

            k[i] = -1.0F * (r[i] + sum) / E[i - 1];             // Equation 38b, Makhoul

            if (Math.abs(k[i]) > 1.0) {
                k[i] = 0.0F;
            }

            a[i][i] = k[i];

            for (j = 1; j <= i - 1; j++) {
                a[i][j] = a[i - 1][j] + k[i] * a[i - 1][i - j];	// Equation 38c, Makhoul

            }

            E[i] = (1.0F - k[i] * k[i]) * E[i - 1];		// Equation 38d, Makhoul
        }

        ak[0] = 1.0F;
        for (i = 1; i <= Defines.LPC_ORD; i++) {
            ak[i] = a[Defines.LPC_ORD][i];
        }

        return ak;
    }
}
