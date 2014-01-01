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
 * Codec2 a speech vocoder. Compresses PCM encoded audio derived from a computer
 * soundcard, and outputs data bits in a packed byte format.
 *
 * <p>
 * The class is designed to work full-duplex, with independent threads each
 * running encode() and decode() functions. Duplicate FFT's and Frame Models are
 * provided for each.
 *
 * <p>
 * NOTE: The audio samples must be sign + 15 bits maximum in an short integer
 * array, with a stable 8 kHz sample rate.
 *
 * <p>
 * Copyright (C) 1990-2013 David Rowe<br>
 * All Rights Reserved
 */
public final class Codec2 {

    public static final int MODE_3200 = 1;
    public static final int MODE_2300 = 2;
    public static final int MODE_1600 = 3;
    public static final int MODE_1400 = 4;
    public static final int MODE_1300 = 5;
    public static final int MODE_1175 = 6;
    public static final int MODES = 6;
    //
    public static final float LPCPF_GAMMA = 0.5F;
    public static final float LPCPF_BETA = 0.2F;
    //
    private static final int TW = 40;                           // Trapezoidal synthesis window overlap
    private static final int LSP_PRED_VQ_INDEXES = 3;
    private static final int SCALAR_INDEXES = 10;
    //
    private final Model[] encodeModels = new Model[4];          // models for encode process
    private final Model[] decodeModels = new Model[4];          // models for decode process
    private Model decodePrev_model;
    //
    private int mode;
    private final Complex[] W;                                  // DFT of w[]
    private Complex[] Sw;                                       // the speech spectrum
    private final float[] w;                                    // time domain hamming window
    private final float[] Sn;                                   // input speech
    private final float[] Sn_;                                  // synthesised output speech
    private final float[] Pn;                                   // trapezoidal synthesis window
    private final float[] prev_lsps_dec;                        // previous frame's LSPs
    private final float[] xq_enc;                               // joint pitch and energy VQ states
    private final float[] xq_dec;
    private final float[] xq_err;                               // running joint pitch and energy error
    private final float[] ex_phase;                             // excitation model phase track
    private final float[] bg_est;                               // background noise estimate for post filter
    private final float[] snr;                                  // running SNR value
    private float prev_Wo_enc;                                  // previous frame's pitch estimate
    private float prev_e_dec;                                   // previous frame's LPC energy

    //
    private float beta;                                         // LPC post filter parameters
    private float gamma;
    //
    private boolean lpc_pf;                                     // LPC post filter enable on/off
    private boolean bass_boost;                                 // LPC post filter bass boost on/off
    /*
     * Classes required to function
     */
    private final FFT fftEncode;
    private final FFT fftDecode;
    private final Nlp nlp;
    private final Sine sine;
    private final Quantise quantise;
    private final Lsp lspc;
    private final Pack pack;

    /*
     * Codec
     * 
     * An instance of this class provides for both encoding and decoding
     */
    public Codec2() {
        int i;
        new Codebookcb();
        new Codebookd();
        new Codebookge();
        new Codebookjvm();

        this.fftEncode = new FFT(Defines.FFT_SIZE);     // used for decoding only
        this.fftDecode = new FFT(Defines.FFT_SIZE);     // used for encoding only
        this.lspc = new Lsp();
        this.quantise = new Quantise(fftDecode, lspc);
        this.sine = new Sine(fftEncode, fftDecode);
        this.nlp = new Nlp(fftEncode);
        this.pack = new Pack();

        this.mode = MODE_1300;

        /*
         * Pre-allocate the four 10 ms voice models
         * for encode and decode (full duplex)
         */
        for (i = 0; i < 4; i++) {
            encodeModels[i] = new Model(this.fftEncode);
            decodeModels[i] = new Model(this.fftDecode);
        }

        this.decodePrev_model = new Model(this.fftDecode);

        this.xq_enc = new float[Codebookge.lsp_cbge[0].k];// used in encode only
        this.xq_dec = new float[Codebookge.lsp_cbge[0].k];// used in decode only
        this.xq_err = new float[Codebookge.lsp_cbge[0].k];// set to 0.0 (size = 2)
        this.bg_est = new float[1];                       // java pointer kludge
        this.ex_phase = new float[1];                     // java pointer kludge
        this.snr = new float[1];                          // ditto
        this.Sn = new float[Defines.M];                   // M = 320 used only in encode
        this.Sn_ = new float[2 * Defines.N];              // N = 80 used only in decode
        this.prev_lsps_dec = new float[Defines.LPC_ORD];  // used in decode only
        this.Pn = new float[2 * Defines.N];               // used in decode only
        this.W = new Complex[Defines.FFT_SIZE];
        this.w = new float[Defines.M];                    // used in encode only

        for (i = 0; i < Defines.FFT_SIZE; i++) {
            this.W[i] = new Complex();
        }

        /*
         * Only need to initialize these once
         */
        make_analysis_window();                           // used during encode
        make_synthesis_window();
    }

    /**
     * Return the last SNR measurement
     *
     * @return a float representing the signal to noise ratio
     */
    public float codec2_getSNR() {
        return this.snr[0];
    }

    /**
     * Return the error value of the WOE method which is used in modes 2300,
     * 1400, and 1175.
     *
     * @return a float array representing the joint VQ pitch and energy errors
     */
    public float[] codec2_getJointVQError() {
        return this.xq_err;
    }

    /**
     * Enable/Disable post filter. Enabled with bass_boost by default in init()
     *
     * @param enable
     * @param bass
     * @param b
     * @param g
     */
    public void codec2_setPostFilter(boolean enable, boolean bass, float b, float g) {
        this.lpc_pf = enable;
        this.bass_boost = bass;
        this.beta = b;
        this.gamma = g;
    }

    /**
     * The codec provides for multiple modes:
     *
     * <p>
     * mode 1 = 3200<br>
     * mode 2 = 2300<br>
     * mode 3 = 1600<br>
     * mode 4 = 1400<br>
     * mode 5 = 1300<br>
     * mode 6 = 1175<br>
     *
     * @param newmode an int selects the codec mode of operation
     * @return an int 0 on success, -1 on failure
     */
    public int codec2_setMode(int newmode) {
        if (newmode > 0 && newmode <= MODES) {
            this.mode = newmode;
            codec2_reset();
            return 0;
        } else {
            return -1;
        }
    }

    /**
     * Returns the number of bits per frame.
     *
     * @return
     */
    public int codec2_getBitsPerFrame() {
        switch (this.mode) {
            case MODE_3200:
            case MODE_1600:
                return 64;
            case MODE_2300:
                return 46;
            case MODE_1175:
                return 47;
            case MODE_1400:
                return 56;
            case MODE_1300:
                return 52;
            default:
                return -1;
        }
    }

    /**
     * Returns the number of bytes per frame.
     *
     * @return
     */
    public int codec2_getBytesPerFrame() {
        return (codec2_getBitsPerFrame() + 7) / 8;
    }

    /**
     * Returns the number of PCM samples per frame.
     *
     * @return an int representing the number of audio samples per frame
     */
    public int codec2_getSamplesPerFrame() {
        switch (this.mode) {
            case MODE_3200:
            case MODE_2300:
                return 160;
            case MODE_1600:
            case MODE_1400:
            case MODE_1300:
            case MODE_1175:
                return 320;
            default:
                return -1;
        }
    }

    /**
     * The output of packed bits from the input speech buffer
     *
     * <p>
     * The byte array is (codec2_getBitsPerFrame() + 7) / 8 bytes long
     *
     * @param bits
     * @param speech
     */
    public void codec2_encode(byte[] bits, short[] speech) {
        /*
         * This assumes the referenced byte array is not initialized before call
         */

        for (int i = 0; i < bits.length; i++) {
            bits[i] = 0;
        }

        switch (this.mode) {
            case MODE_3200:
                encode_3200(bits, speech);
                break;
            case MODE_2300:
                encode_2300(bits, speech);
                break;
            case MODE_1600:
                encode_1600(bits, speech);
                break;
            case MODE_1400:
                encode_1400(bits, speech);
                break;
            case MODE_1300:
                encode_1300(bits, speech);
                break;
            case MODE_1175:
                encode_1175(bits, speech);
        }
    }

    /**
     * The output of a speech buffer from an input of packed bits
     *
     * @param speech
     * @param bits
     */
    public void codec2_decode(short[] speech, byte[] bits) {
        switch (this.mode) {
            case MODE_3200:
                decode_3200(speech, bits);
                break;
            case MODE_2300:
                decode_2300(speech, bits);
                break;
            case MODE_1600:
                decode_1600(speech, bits);
                break;
            case MODE_1400:
                decode_1400(speech, bits);
                break;
            case MODE_1300:
                decode_1300(speech, bits);
                break;
            case MODE_1175:
                decode_1175(speech, bits);
        }
    }

    /**
     * This provides the user a way to reset everything back to defaults in the
     * current mode.
     */
    private void codec2_reset() {
        int i;

        this.decodePrev_model.reInit();
        this.decodePrev_model.setWo((Defines.TWO_PI / Defines.P_MAX));    // 0.039269875
        this.decodePrev_model.setL((int) (Defines.TWO_PI / Defines.P_MAX));       // 0 ??

        for (i = 0; i < 4; i++) {
            decodeModels[i].reInit();
            encodeModels[i].reInit();
        }

        this.bg_est[0] = 0.0F;                            // used in decode only
        this.ex_phase[0] = 0.0F;                          // used in decode only
        this.snr[0] = 0.0F;                               // Signal/Noise Ratio

        this.prev_e_dec = 1.0F;                           // used in decode only
        this.prev_Wo_enc = 0.0F;                          // used in encode only

        for (i = 0; i < Defines.M; i++) {
            this.Sn[i] = 1.0F;
        }

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = i * ((float) Math.PI / (Defines.LPC_ORD + 1));
        }

        codec2_setPostFilter(true, true, Codec2.LPCPF_BETA, Codec2.LPCPF_GAMMA);
    }

    /*---------------------------------------------------------------------------*\

     Encodes 160 speech samples (20ms of speech) into 64 bits.

     The codec2 algorithm actually operates internally on 10ms (80
     sample) frames, so we run the encoding algorithm twice.  On the
     first frame we just send the voicing bits.  On the second frame we
     send all model parameters.  Compared to 2400 we use a larger number
     of bits for the LSPs and non-VQ pitch and energy.

     The bit allocation is:

     Parameter                      bits/frame
     --------------------------------------
     Harmonic magnitudes (LSPs)     50
     Pitch (Wo)                      7
     Energy                          5
     Voicing (10ms update)           2
     TOTAL                          64

     \*---------------------------------------------------------------------------*/
    private void encode_3200(byte[] bits, short[] speech) {
        float[] lsps = new float[Defines.LPC_ORD];
        int[] lspd_indexes = new int[Defines.LPC_ORD];

        pack.Init();
        encodeModels[0].reInit();

        /* first 10ms analysis frame - we just want voicing */
        analyse_one_frame(encodeModels[0], speech, 0);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* second 10ms analysis frame */
        analyse_one_frame(encodeModels[0], speech, Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        pack.pack(bits, quantise.encode_Wo(encodeModels[0].getWo()), Defines.WO_BITS);
        pack.pack(bits, quantise.encode_energy(quantise.speech_to_uq_lsps(lsps, this.Sn, this.w)), Defines.E_BITS);

        quantise.encode_lspds_scalar(lspd_indexes, lsps);

        for (int i = 0; i < SCALAR_INDEXES; i++) {
            pack.pack(bits, lspd_indexes[i], quantise.lspd_bits(i));
        }
    }

    /*---------------------------------------------------------------------------*\

     Decodes a frame of 64 bits into 160 samples (20ms) of speech.

     \*---------------------------------------------------------------------------*/
    private void decode_3200(short[] speech, byte[] bits) {
        float[][] lsps = new float[2][Defines.LPC_ORD];
        float[][] ak = new float[2][Defines.LPC_ORD + 1];
        float[] e = new float[2];
        int[] lspd_indexes = new int[Defines.LPC_ORD];
        int i, Wo_index, e_index;

        pack.Init();

        for (i = 0; i < 2; i++) {
            decodeModels[i].reInit();
        }

        /* unpack bits from channel ------------------------------------*/

        /*
         * this will partially fill the model params for the 2 x 10ms
         * frames
         */
        decodeModels[0].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[1].setVoiced((pack.unpack(bits, 1) == 1));

        Wo_index = pack.unpack(bits, Defines.WO_BITS);
        decodeModels[1].setWo(quantise.decode_Wo(Wo_index));
        decodeModels[1].setL((int) (Math.PI / decodeModels[1].getWo()));

        e_index = pack.unpack(bits, Defines.E_BITS);
        e[1] = quantise.decode_energy(e_index);

        for (i = 0; i < SCALAR_INDEXES; i++) {
            lspd_indexes[i] = pack.unpack(bits, quantise.lspd_bits(i));
        }

        quantise.decode_lspds_scalar(lsps[1], lspd_indexes);

        /*
         * Wo and energy are sampled every 20ms, so we interpolate just 1
         * 10ms frame between 20ms samples
         */
        interp_Wo(decodeModels[0], decodePrev_model, decodeModels[1]);      // interpolate
        e[0] = interp_energy(this.prev_e_dec, e[1]);

        /*
         * LSPs are sampled every 20ms so we interpolate the frame in
         * between, then recover spectral amplitudes
         */
        interpolate_lsp(lsps[0], this.prev_lsps_dec, lsps[1], 0.5F);

        for (i = 0; i < 2; i++) {
            lspc.lsp_to_lpc(lsps[i], ak[i]);
            quantise.aks_to_M2(ak[i], decodeModels[i], e[i], this.lpc_pf, this.bass_boost,
                    this.beta, this.gamma, this.snr);
            decodeModels[i].apply_lpc_correction();
            synthesise_one_frame(decodeModels[i], speech, Defines.N * i, ak[i]);      // synthesise
        }

        // update memories for next frame
        this.decodePrev_model = decodeModels[1];
        this.prev_e_dec = e[1];

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = lsps[1][i];
        }
    }

    /*---------------------------------------------------------------------------*\

     Encodes 160 speech samples (20ms of speech) into 46 bits.

     The codec2 algorithm actually operates internally on 10ms (80
     sample) frames, so we run the encoding algorithm twice.  On the
     first frame we just send the voicing bit.  On the second frame we
     send all model parameters.

     The bit allocation is:

     Parameter                      bits/frame
     --------------------------------------
     Harmonic magnitudes (LSPs)     36
     Joint VQ of Energy and Wo       8
     Voicing (10ms update)           2
     TOTAL                          46

     \*---------------------------------------------------------------------------*/
    private void encode_2300(byte[] bits, short[] speech) {
        float[] lsps = new float[Defines.LPC_ORD];
        int[] lsp_indexes = new int[Defines.LPC_ORD];

        pack.Init();

        encodeModels[0].reInit();

        analyse_one_frame(encodeModels[0], speech, 0);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* second 10ms analysis frame */
        analyse_one_frame(encodeModels[0], speech, Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);
        pack.pack(bits, quantise.encode_WoE(encodeModels[0],
                quantise.speech_to_uq_lsps(lsps, this.Sn, this.w), this.xq_err, this.xq_enc), Defines.WO_E_BITS);

        quantise.encode_lsps_scalar(lsp_indexes, lsps);

        for (int i = 0; i < SCALAR_INDEXES; i++) {
            pack.pack(bits, lsp_indexes[i], quantise.lsp_bits(i));
        }
    }

    /*---------------------------------------------------------------------------*\

     Decodes frames of 46 bits into 160 samples (20ms) of speech.

     \*---------------------------------------------------------------------------*/
    private void decode_2300(short[] speech, byte[] bits) {
        int[] lsp_indexes = new int[Defines.LPC_ORD];
        float[][] lsps = new float[2][Defines.LPC_ORD];
        int i, WoE_index;
        float[] e = new float[2];
        float[][] ak = new float[2][Defines.LPC_ORD + 1];

        pack.Init();

        for (i = 0; i < 2; i++) {
            decodeModels[i].reInit();
        }

        /* unpack bits from channel ------------------------------------*/

        /*
         * this will partially fill the model params for the 2 x 10ms
         * frames
         */
        decodeModels[0].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[1].setVoiced((pack.unpack(bits, 1) == 1));

        WoE_index = pack.unpack(bits, Defines.WO_E_BITS);
        e[1] = quantise.decode_WoE(decodeModels[1], WoE_index, this.xq_dec);

        for (i = 0; i < SCALAR_INDEXES; i++) {
            lsp_indexes[i] = pack.unpack(bits, quantise.lsp_bits(i));
        }

        quantise.decode_lsps_scalar(lsps[1], lsp_indexes);

        /* interpolate ------------------------------------------------*/

        /*
         * Wo and energy are sampled every 20ms, so we interpolate just one
         * 10ms frame between 20ms samples
         */
        interp_Wo(decodeModels[0], decodePrev_model, decodeModels[1]);
        e[0] = interp_energy(this.prev_e_dec, e[1]);

        /*
         * LSPs are sampled every 20ms so we interpolate the frame in
         * between, then recover spectral amplitudes
         */
        interpolate_lsp(lsps[0], this.prev_lsps_dec, lsps[1], 0.5F);

        for (i = 0; i < 2; i++) {
            lspc.lsp_to_lpc(lsps[i], ak[i]);
            quantise.aks_to_M2(ak[i], decodeModels[i], e[i], this.lpc_pf, this.bass_boost,
                    this.beta, this.gamma, this.snr);
            decodeModels[i].apply_lpc_correction();
            synthesise_one_frame(decodeModels[i], speech, Defines.N * i, ak[i]);       // synthesise
        }

        /* update memories for next frame ----------------------------*/
        this.decodePrev_model = decodeModels[1];
        this.prev_e_dec = e[1];

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = lsps[1][i];
        }
    }

    /*---------------------------------------------------------------------------*\

     Encodes 320 speech samples (40ms of speech) into 64 bits.

     The codec2 algorithm actually operates internally on 10ms (80
     sample) frames, so we run the encoding algorithm 4 times:

     frame 1: voicing bit
     frame 2: voicing bit, Wo and E
     frame 3: voicing bit
     frame 4: voicing bit, Wo and E, scalar LSPs

     The bit allocation is:

     Parameter                      frame 2  frame 4   Total
     -------------------------------------------------------
     Harmonic magnitudes (LSPs)      0       36        36
     Pitch (Wo)                      7        7        14
     Energy                          5        5        10
     Voicing (10ms update)           2        2         4
     TOTAL                          14       50        64

     \*---------------------------------------------------------------------------*/
    private void encode_1600(byte[] bits, short[] speech) {
        float[] lsps = new float[Defines.LPC_ORD];
        int[] lsp_indexes = new int[Defines.LPC_ORD];

        pack.Init();

        encodeModels[0].reInit();

        /*
         * frame 1: - voicing
         */
        analyse_one_frame(encodeModels[0], speech, 0);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);     // 1 bit

        /*
         * frame 2: - voicing, scalar Wo & E
         */
        analyse_one_frame(encodeModels[0], speech, Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);                    // 1 bit
        pack.pack(bits, quantise.encode_Wo(encodeModels[0].getWo()), Defines.WO_BITS);  // 7 bits

        /* need to run this just to get LPC energy */
        pack.pack(bits, quantise.encode_energy(quantise.speech_to_uq_lsps(lsps, this.Sn, this.w)), Defines.E_BITS);         // 5 bits

        /*
         * frame 3: - voicing
         */
        analyse_one_frame(encodeModels[0], speech, 2 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);     // 1 bit

        /*
         * frame 4: - voicing, scalar Wo & E, scalar LSPs
         */
        analyse_one_frame(encodeModels[0], speech, 3 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);     // 1 bit
        pack.pack(bits, quantise.encode_Wo(encodeModels[0].getWo()), Defines.WO_BITS);       // 7 bits
        pack.pack(bits, quantise.encode_energy(quantise.speech_to_uq_lsps(lsps, this.Sn, this.w)), Defines.E_BITS); // 5 bits

        quantise.encode_lsps_scalar(lsp_indexes, lsps);

        /*
         * Codebook has 10 entries, with 36 bits
         */
        for (int i = 0; i < SCALAR_INDEXES; i++) {      // 36 bits
            pack.pack(bits, lsp_indexes[i], quantise.lsp_bits(i));
        }
    }

    /*---------------------------------------------------------------------------*\

     Decodes frames of 64 bits into 320 samples (40ms) of speech.

     \*---------------------------------------------------------------------------*/
    private void decode_1600(short[] speech, byte[] bits) {
        int[] lsp_indexes = new int[Defines.LPC_ORD];
        float[][] lsps = new float[4][Defines.LPC_ORD];
        int i, Wo_index, e_index;
        float[] e = new float[4];
        float[][] ak = new float[4][Defines.LPC_ORD + 1];
        float weight;

        pack.Init();

        for (i = 0; i < 4; i++) {
            decodeModels[i].reInit();
        }

        // unpack bits from channel
        // this will partially fill the model params for the 4 x 10ms frames
        decodeModels[0].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[1].setVoiced((pack.unpack(bits, 1) == 1));

        Wo_index = pack.unpack(bits, Defines.WO_BITS);
        decodeModels[1].setWo(quantise.decode_Wo(Wo_index));
        decodeModels[1].setL((int) (Math.PI / decodeModels[1].getWo()));

        e_index = pack.unpack(bits, Defines.E_BITS);
        e[1] = quantise.decode_energy(e_index);

        decodeModels[2].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[3].setVoiced((pack.unpack(bits, 1) == 1));

        Wo_index = pack.unpack(bits, Defines.WO_BITS);
        decodeModels[3].setWo(quantise.decode_Wo(Wo_index));
        decodeModels[3].setL((int) (Math.PI / decodeModels[3].getWo()));

        e_index = pack.unpack(bits, Defines.E_BITS);
        e[3] = quantise.decode_energy(e_index);

        for (i = 0; i < SCALAR_INDEXES; i++) {
            lsp_indexes[i] = pack.unpack(bits, quantise.lsp_bits(i));
        }

        quantise.decode_lsps_scalar(lsps[3], lsp_indexes);

        /* interpolate ------------------------------------------------*/

        /*
         * Wo and energy are sampled every 20ms, so we interpolate just one
         * 10ms frame between 20ms samples
         */
        interp_Wo(decodeModels[0], decodePrev_model, decodeModels[1]);
        e[0] = interp_energy(this.prev_e_dec, e[1]);
        interp_Wo(decodeModels[2], decodeModels[1], decodeModels[3]);
        e[2] = interp_energy(e[1], e[3]);

        /*
         * LSPs are sampled every 40ms so we interpolate the 3 frames in
         * between, then recover spectral amplitudes
         */
        for (i = 0, weight = 0.25F; i < 3; i++, weight += 0.25F) {
            interpolate_lsp(lsps[i], this.prev_lsps_dec, lsps[3], weight);
        }

        for (i = 0; i < 4; i++) {
            lspc.lsp_to_lpc(lsps[i], ak[i]);
            quantise.aks_to_M2(ak[i], decodeModels[i], e[i], this.lpc_pf, this.bass_boost,
                    this.beta, this.gamma, this.snr);
            decodeModels[i].apply_lpc_correction();
            synthesise_one_frame(decodeModels[i], speech, Defines.N * i, ak[i]);
        }

        /* update memories for next frame ----------------------------*/
        this.decodePrev_model = decodeModels[3];
        this.prev_e_dec = e[3];

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = lsps[3][i];
        }
    }

    /*---------------------------------------------------------------------------*\

     Encodes 320 speech samples (40ms of speech) into 56 bits.

     The codec2 algorithm actually operates internally on 10ms (80
     sample) frames, so we run the encoding algorithm 4 times:

     frame 0: voicing bit
     frame 1: voicing bit, joint VQ of Wo and E
     frame 2: voicing bit
     frame 3: voicing bit, joint VQ of Wo and E, scalar LSPs

     The bit allocation is:

     Parameter                      frame 2  frame 4   Total
     -------------------------------------------------------
     Harmonic magnitudes (LSPs)      0       36        36
     Energy+Wo                       8        8        16
     Voicing (10ms update)           2        2         4
     TOTAL                          10       46        56

     \*---------------------------------------------------------------------------*/
    private void encode_1400(byte[] bits, short[] speech) {
        float[] lsps = new float[Defines.LPC_ORD];
        int[] lsp_indexes = new int[Defines.LPC_ORD];

        pack.Init();

        encodeModels[0].reInit();

        /* frame 1: - voicing ---------------------------------------------*/
        analyse_one_frame(encodeModels[0], speech, 0);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* frame 2: - voicing, joint Wo & E -------------------------------*/
        analyse_one_frame(encodeModels[0], speech, Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);
        pack.pack(bits, quantise.encode_WoE(encodeModels[0],
                quantise.speech_to_uq_lsps(lsps, this.Sn, this.w), this.xq_err, this.xq_enc), Defines.WO_E_BITS);

        /* frame 3: - voicing ---------------------------------------------*/
        analyse_one_frame(encodeModels[0], speech, 2 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* frame 4: - voicing, joint Wo & E, scalar LSPs ------------------*/
        analyse_one_frame(encodeModels[0], speech, 3 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);
        pack.pack(bits, quantise.encode_WoE(encodeModels[0],
                quantise.speech_to_uq_lsps(lsps, this.Sn, this.w), this.xq_err, this.xq_enc), Defines.WO_E_BITS);

        quantise.encode_lsps_scalar(lsp_indexes, lsps);

        for (int i = 0; i < SCALAR_INDEXES; i++) {
            pack.pack(bits, lsp_indexes[i], quantise.lsp_bits(i));
        }
    }

    /*---------------------------------------------------------------------------*\

     Decodes frames of 56 bits into 320 samples (40ms) of speech.

     \*---------------------------------------------------------------------------*/
    private void decode_1400(short[] speech, byte[] bits) {
        int[] lsp_indexes = new int[Defines.LPC_ORD];
        float[][] lsps = new float[4][Defines.LPC_ORD];
        int i, WoE_index;
        float[] e = new float[4];
        float[][] ak = new float[4][Defines.LPC_ORD + 1];
        float weight;

        pack.Init();

        for (i = 0; i < 4; i++) {
            decodeModels[i].reInit();
        }

        /* unpack bits from channel ------------------------------------*/

        /*
         * this will partially fill the model params for the 4 x 10ms
         * frames
         */
        decodeModels[0].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[1].setVoiced((pack.unpack(bits, 1) == 1));

        WoE_index = pack.unpack(bits, Defines.WO_E_BITS);
        e[1] = quantise.decode_WoE(decodeModels[1], WoE_index, this.xq_dec);

        decodeModels[2].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[3].setVoiced((pack.unpack(bits, 1) == 1));

        WoE_index = pack.unpack(bits, Defines.WO_E_BITS);
        e[3] = quantise.decode_WoE(decodeModels[3], WoE_index, this.xq_dec);

        for (i = 0; i < SCALAR_INDEXES; i++) {
            lsp_indexes[i] = pack.unpack(bits, quantise.lsp_bits(i));
        }

        quantise.decode_lsps_scalar(lsps[3], lsp_indexes);

        /* interpolate ------------------------------------------------*/

        /*
         * Wo and energy are sampled every 20ms, so we interpolate just 1
         * 10ms frame between 20ms samples
         */
        interp_Wo(decodeModels[0], decodePrev_model, decodeModels[1]);
        e[0] = interp_energy(this.prev_e_dec, e[1]);

        interp_Wo(decodeModels[2], decodeModels[1], decodeModels[3]);
        e[2] = interp_energy(e[1], e[3]);

        /*
         * LSPs are sampled every 40ms so we interpolate the 3 frames in
         * between, then recover spectral amplitudes
         */
        for (i = 0, weight = 0.25F; i < 3; i++, weight += 0.25F) {
            interpolate_lsp(lsps[i], this.prev_lsps_dec, lsps[3], weight);
        }

        for (i = 0; i < 4; i++) {
            lspc.lsp_to_lpc(lsps[i], ak[i]);
            quantise.aks_to_M2(ak[i], decodeModels[i], e[i], this.lpc_pf, this.bass_boost,
                    this.beta, this.gamma, this.snr);
            decodeModels[i].apply_lpc_correction();
            synthesise_one_frame(decodeModels[i], speech, Defines.N * i, ak[i]);
        }

        /* update memories for next frame ----------------------------*/
        this.decodePrev_model = decodeModels[3];
        this.prev_e_dec = e[3];

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = lsps[3][i];
        }
    }

    /*---------------------------------------------------------------------------*\

     Encodes 320 speech samples (40ms of speech) into 52 bits.

     The codec2 algorithm actually operates internally on 10ms (80
     sample) frames, so we run the encoding algorithm 4 times:

     frame 0: voicing bit
     frame 1: voicing bit,
     frame 2: voicing bit
     frame 3: voicing bit, Wo and E, scalar LSPs

     The bit allocation is:

     Parameter                      frame 2  frame 4   Total
     -------------------------------------------------------
     Harmonic magnitudes (LSPs)      0       36        36
     Pitch (Wo)                      0        7         7
     Energy                          0        5         5
     Voicing (10ms update)           2        2         4
     TOTAL                           2       50        52

     \*---------------------------------------------------------------------------*/
    private void encode_1300(byte[] bits, short[] speech) {
        float[] lsps = new float[Defines.LPC_ORD];
        int[] lsp_indexes = new int[Defines.LPC_ORD];

        pack.Init();

        encodeModels[0].reInit();

        /* first 10ms analysis frame - we just want voicing */
        analyse_one_frame(encodeModels[0], speech, 0);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* second 10ms analysis frame */
        analyse_one_frame(encodeModels[0], speech, Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* frame 3: - voicing ---------------------------------------------*/
        analyse_one_frame(encodeModels[0], speech, 2 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* frame 4: - voicing, scalar Wo & E, scalar LSPs ------------------*/
        analyse_one_frame(encodeModels[0], speech, 3 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);
        pack.pack(bits, quantise.encode_Wo(encodeModels[0].getWo()), Defines.WO_BITS);  // 7 bits
        pack.pack(bits, quantise.encode_energy(quantise.speech_to_uq_lsps(lsps, this.Sn, this.w)), Defines.E_BITS); // 5 bits

        quantise.encode_lsps_scalar(lsp_indexes, lsps);

        for (int i = 0; i < SCALAR_INDEXES; i++) {
            pack.pack(bits, lsp_indexes[i], quantise.lsp_bits(i));      // 36 bits total
        }
    }

    /*---------------------------------------------------------------------------*\

     Decodes frames of 52 bits into 320 samples (40ms) of speech.

     \*---------------------------------------------------------------------------*/
    private void decode_1300(short[] speech, byte[] bits) {
        int[] lsp_indexes = new int[Defines.LPC_ORD];
        float[][] lsps = new float[4][Defines.LPC_ORD];
        int i, Wo_index, e_index;
        float[] e = new float[4];
        float[][] ak = new float[4][Defines.LPC_ORD + 1];
        float weight;

        pack.Init();

        for (i = 0; i < 4; i++) {
            decodeModels[i].reInit();
        }

        /* unpack bits from channel ------------------------------------*/

        /*
         * this will partially fill the model params for the 4 x 10ms
         * frames
         */
        decodeModels[0].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[1].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[2].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[3].setVoiced((pack.unpack(bits, 1) == 1));

        Wo_index = pack.unpack(bits, Defines.WO_BITS);
        decodeModels[3].setWo(quantise.decode_Wo(Wo_index));
        decodeModels[3].setL((int) (Math.PI / decodeModels[3].getWo()));
        e_index = pack.unpack(bits, Defines.E_BITS);
        e[3] = quantise.decode_energy(e_index);

        for (i = 0; i < SCALAR_INDEXES; i++) {
            lsp_indexes[i] = pack.unpack(bits, quantise.lsp_bits(i));
        }

        quantise.decode_lsps_scalar(lsps[3], lsp_indexes);

        /* interpolate ------------------------------------------------*/

        /*
         * Wo, energy, and LSPs are sampled every 40ms so we interpolate
         * the 3 frames in between
         */
        for (i = 0, weight = 0.25F; i < 3; i++, weight += 0.25F) {
            interpolate_lsp(lsps[i], this.prev_lsps_dec, lsps[3], weight);
            interp_Wo2(decodeModels[i], decodePrev_model, decodeModels[3], weight);
            e[i] = interp_energy2(this.prev_e_dec, e[3], weight);
        }

        /*
         * then recover spectral amplitudes
         */
        for (i = 0; i < 4; i++) {
            lspc.lsp_to_lpc(lsps[i], ak[i]);
            quantise.aks_to_M2(ak[i], decodeModels[i], e[i], this.lpc_pf, this.bass_boost,
                    this.beta, this.gamma, this.snr);
            decodeModels[i].apply_lpc_correction();
            synthesise_one_frame(decodeModels[i], speech, Defines.N * i, ak[i]);
        }

        /* update memories for next frame ----------------------------*/
        this.decodePrev_model = decodeModels[3];
        this.prev_e_dec = e[3];

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = lsps[3][i];
        }
    }

    /*---------------------------------------------------------------------------*\

     Encodes 320 speech samples (40ms of speech) into 47 bits.  

     The codec2 algorithm actually operates internally on 10ms (80
     sample) frames, so we run the encoding algorithm four times:

     frame 0: voicing bit
     frame 1: voicing bit, joint VQ of Wo and E
     frame 2: voicing bit
     frame 3: voicing bit, joint VQ of Wo and E, VQ LSPs

     The bit allocation is:

     Parameter                      frame 2  frame 4   Total
     -------------------------------------------------------
     Harmonic magnitudes (LSPs)      0       27        27
     Energy+Wo                       8        8        16
     Voicing (10ms update)           2        2         4
     TOTAL                          10       37        47
 
     \*---------------------------------------------------------------------------*/
    private void encode_1175(byte[] bits, short[] speech) {
        float[] lsps = new float[Defines.LPC_ORD];
        int[] lsp_indexes = new int[Defines.LPC_ORD];

        pack.Init();

        /* first 10ms analysis frame - we just want voicing */
        analyse_one_frame(encodeModels[0], speech, 0);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* second 10ms analysis frame */
        analyse_one_frame(encodeModels[0], speech, Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* need to run this just to get LPC energy */
        pack.pack(bits, quantise.encode_WoE(encodeModels[0],
                quantise.speech_to_uq_lsps(lsps, this.Sn, this.w), this.xq_err, this.xq_enc), Defines.WO_E_BITS);

        /* frame 3: - voicing ---------------------------------------------*/
        analyse_one_frame(encodeModels[0], speech, 2 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);

        /* frame 4: - voicing, joint Wo & E, scalar LSPs ------------------*/
        analyse_one_frame(encodeModels[0], speech, 3 * Defines.N);
        pack.pack(bits, ((encodeModels[0].getVoiced()) ? 1 : 0), 1);
        pack.pack(bits, quantise.encode_WoE(encodeModels[0],
                quantise.speech_to_uq_lsps(lsps, this.Sn, this.w), this.xq_err, this.xq_enc), Defines.WO_E_BITS);

        quantise.encode_lsps_vq(lsp_indexes, lsps);

        for (int i = 0; i < LSP_PRED_VQ_INDEXES; i++) {
            pack.pack(bits, lsp_indexes[i], quantise.lsp_pred_vq_bits(i));
        }
    }


    /*---------------------------------------------------------------------------*\

     Decodes frames of 47 bits into 320 samples (40ms) of speech.

     \*---------------------------------------------------------------------------*/
    private void decode_1175(short[] speech, byte[] bits) {
        int[] lsp_indexes = new int[Defines.LPC_ORD];
        float[][] lsps = new float[4][Defines.LPC_ORD];
        int i, WoE_index;
        float[] e = new float[4];
        float[][] ak = new float[4][Defines.LPC_ORD + 1];
        float weight;

        pack.Init();

        for (i = 0; i < 4; i++) {
            decodeModels[i].reInit();
        }

        /*
         * unpack bits from channel
         */
        decodeModels[0].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[1].setVoiced((pack.unpack(bits, 1) == 1));

        WoE_index = pack.unpack(bits, Defines.WO_E_BITS);
        e[1] = quantise.decode_WoE(decodeModels[1], WoE_index, this.xq_dec);

        decodeModels[2].setVoiced((pack.unpack(bits, 1) == 1));
        decodeModels[3].setVoiced((pack.unpack(bits, 1) == 1));

        WoE_index = pack.unpack(bits, Defines.WO_E_BITS);
        e[3] = quantise.decode_WoE(decodeModels[3], WoE_index, this.xq_dec);

        for (i = 0; i < LSP_PRED_VQ_INDEXES; i++) {
            lsp_indexes[i] = pack.unpack(bits, quantise.lsp_pred_vq_bits(i));
        }

        quantise.decode_lsps_vq(lsps[3], lsp_indexes);

        /* interpolate ------------------------------------------------*/

        /*
         * Wo and energy are sampled every 20ms, so we interpolate
         * just one 10ms frame between 20ms samples
         */
        interp_Wo(decodeModels[0], decodePrev_model, decodeModels[1]);
        e[0] = interp_energy(this.prev_e_dec, e[1]);

        interp_Wo(decodeModels[2], decodeModels[1], decodeModels[3]);
        e[2] = interp_energy(e[1], e[3]);

        /*
         * LSPs are sampled every 40ms so we interpolate the 3 frames in
         * between, then recover spectral amplitudes
         */
        for (i = 0, weight = 0.25F; i < 3; i++, weight += 0.25F) {
            interpolate_lsp(lsps[i], this.prev_lsps_dec, lsps[3], weight);
        }

        for (i = 0; i < 4; i++) {
            lspc.lsp_to_lpc(lsps[i], ak[i]);
            quantise.aks_to_M2(ak[i], decodeModels[i], e[i], this.lpc_pf, this.bass_boost,
                    this.beta, this.gamma, this.snr);
            decodeModels[i].apply_lpc_correction();
            synthesise_one_frame(decodeModels[i], speech, Defines.N * i, ak[i]);
        }

        /* update memories for next frame ----------------------------*/
        this.decodePrev_model = decodeModels[3];
        this.prev_e_dec = e[3];

        for (i = 0; i < Defines.LPC_ORD; i++) {
            this.prev_lsps_dec[i] = lsps[3][i];
        }
    }

    /**
     * Initialize window synthesis function
     *
     * Generates the trapezoidal (Parzen) synthesis window.
     */
    private void make_synthesis_window() {
        int i;

        /* Generate Parzen window in time domain */
        for (i = 0; i < Defines.N / 2 - Codec2.TW; i++) {
            this.Pn[i] = 0.0F;
        }

        float win = 0.0F;
        for (i = Defines.N / 2 - Codec2.TW; i < Defines.N / 2 + Codec2.TW; win += 1.0 / (2 * Codec2.TW), i++) {
            this.Pn[i] = win;
        }

        for (i = Defines.N / 2 + Codec2.TW; i < 3 * Defines.N / 2 - Codec2.TW; i++) {
            this.Pn[i] = 1.0F;
        }

        win = 1.0F;
        for (i = 3 * Defines.N / 2 - Codec2.TW; i < 3 * Defines.N / 2 + Codec2.TW; win -= 1.0 / (2 * Codec2.TW), i++) {
            this.Pn[i] = win;
        }

        for (i = 3 * Defines.N / 2 + Codec2.TW; i < 2 * Defines.N; i++) {
            this.Pn[i] = 0.0F;
        }
    }

    /**
     * Initialize window analysis function
     *
     * Generates the time domain analysis window and it's DFT.
     */
    private void make_analysis_window() {
        int i, j;

        /*
         * Generate Hamming window centered on M-sample pitch analysis window
         * 
         *         0            M/2           M-1
         *         |-------------|-------------|
         *               |-------|-------|
         *                  NW samples
         *
         * All our analysis/synthsis is centred on the M/2 sample.
         */
        for (i = 0; i < (Defines.M / 2) - (Defines.NW / 2); i++) {
            this.w[i] = 0.0F;
        }

        float m = 0.0F;

        for (i = (Defines.M / 2) - (Defines.NW / 2), j = 0; i < (Defines.M / 2) + (Defines.NW / 2); i++, j++) {
            this.w[i] = (float) (0.5 - 0.5 * Math.cos(Defines.TWO_PI * (double) j / (Defines.NW - 1)));
            m += (this.w[i] * this.w[i]);
        }

        for (i = Defines.M / 2 + Defines.NW / 2; i < Defines.M; i++) {
            this.w[i] = 0.0F;
        }

        /*
         * Normalize - makes freq domain amplitude estimation straight
         * forward
         */
        m = (float) (1.0 / Math.sqrt((double) m * (double) Defines.FFT_SIZE));

        for (i = 0; i < Defines.M; i++) {
            this.w[i] *= m;
        }

        /*
         * Generate DFT of analysis window, used for later processing.  Note
         * we modulo FFT_SIZE shift the time domain window w[], this makes the
         * imaginary part of the DFT W[] equal to zero as the shifted w[] is
         * even about the n=0 time axis if NW is odd.  Having the imag part
         * of the DFT W[] makes computation easier.
         * 
         *          0                      FFT_SIZE-1
         *          |-------------------------|
         * 
         *           ----\               /----
         *                \             /
         *                 \           /      <- shifted version of window w[n]
         *                  \         /
         *                   \       /
         *                    -------
         * 
         *          |---------|     |---------|
         *            NW/2              NW/2
         */
        for (i = 0; i < Defines.FFT_SIZE; i++) {
            this.W[i] = new Complex();
        }

        for (i = 0; i < Defines.NW / 2; i++) {
            this.W[i] = new Complex(this.w[i + Defines.M / 2], 0.0F);
        }

        for (i = Defines.FFT_SIZE - Defines.NW / 2, j = Defines.M / 2 - Defines.NW / 2; i < Defines.FFT_SIZE; i++, j++) {
            this.W[i] = new Complex(this.w[j], 0.0F);
        }

        fftEncode.transform(this.W);

        /*
         * Re-arrange W[] to be symmetrical about FFT_SIZE/2.  Makes later
         * analysis convenient.
         * 
         *         Before:
         * 
         * 
         *          0                 FFT_SIZE-1
         *          |----------|---------|
         *          __                   _
         *            \                 /
         *             \_______________/
         * 
         * After:
         * 
         *          0                 FFT_SIZE-1
         *          |----------|---------|
         *                    ___
         *                   /   \
         *          ________/     \_______
         * 
         */
        Complex temp;

        for (i = 0; i < Defines.FFT_SIZE / 2; i++) {
            temp = this.W[i];
            this.W[i] = this.W[i + Defines.FFT_SIZE / 2];
            this.W[i + Defines.FFT_SIZE / 2] = temp;
        }
    }

    /**
     * Synthesize 80 speech samples (10ms) from model parameters. Limits output
     * level to protect ears when there are bit errors or the input is over
     * driven.
     *
     * This doesn't correct or mask bit errors, just reduces the worst of their
     * damage.
     *
     * @param speech
     * @param index
     * @param model
     * @param ak
     */
    private void synthesise_one_frame(Model model, short[] speech, int index, float[] ak) {
        float gain;
        int i;

        model.phase_synth_zero_order(ak, this.ex_phase);
        model.postfilter(bg_est);
        sine.synthesise(model, this.Sn_, this.Pn);

        /*
         * find maximum sample in frame for ear protection
         */
        float max_sample = 0.0F;

        for (i = 0; i < Defines.N; i++) {
            if (this.Sn_[i] > max_sample) {
                max_sample = this.Sn_[i];
            }
        }

        /*
         * determine how far above set point
         */
        float over = max_sample / 30000.0F;

        /*
         * If we are x dB over set point we reduce level by 2x dB, this
         * attenuates major excursions in amplitude (likely to be caused
         * by bit errors) more than smaller ones
         */
        if (over > 1.0F) {
            gain = 1.0F / (over * over);

            for (i = 0; i < Defines.N; i++) {
                this.Sn_[i] *= gain;
            }
        }

        for (i = 0; i < Defines.N; i++) {
            if (this.Sn_[i] > 32767.0F) {
                speech[i + index] = 32767;
            } else if (this.Sn_[i] < -32767.0F) {
                speech[i + index] = -32767;
            } else {
                speech[i + index] = (short) this.Sn_[i];
            }
        }
    }

    /**
     * Extract sinusoidal model parameters from input speech samples. 10
     * milliseconds of new speech samples are added during each call.
     *
     * @param model
     * @param speech
     * @param index
     */
    private void analyse_one_frame(Model model, short[] speech, int index) {
        float Wo;
        int i;

        /*
         * Prepare for a new 80 samples.
         * 
         * The Sn array is initialized to all 1.0 values in init()
         */
        System.arraycopy(this.Sn, Defines.N, this.Sn, 0, Defines.M - Defines.N);    // M = 320, N = 80, M-N = 240

        /*
         * Now add the new samples to the end as floating point values
         */
        for (i = 0; i < Defines.N; i++) {                   // N = 80
            this.Sn[i + (Defines.M - Defines.N)] = (float) speech[i + index];
        }

        /*
         * Perform an FFT on the Anaysis Window
         */
        sine.dft_speech(this.Sn, this.w);
        this.Sw = sine.getEncodeDomain();                   // results in 512 complex frequency bins

        // Estimate pitch
        nlp.nlp(this.Sn, this.prev_Wo_enc);                 // performs another FFT to determine fundamental
        Wo = (Defines.TWO_PI / nlp.getPitch());
        model.setWo(Wo);
        model.setL((int) Wo);

        // estimate model parameters
        sine.two_stage_pitch_refinement(model, this.Sw);     // refine the pitch
        sine.estimate_amplitudes(model, this.Sw, true);     // estimate phase (false for embedded cpu)
        sine.est_voicing_mbe(model, this.Sw, this.W);

        this.prev_Wo_enc = model.getWo();     // used in NLP post_process_sub_multiples()
    }

    /**
     * Interpolates center 10ms sample of Wo and L samples given two samples
     * 20ms apart. Assumes voicing is available for center (interpolated) frame.
     *
     * @param interp
     * @param prev
     * @param next
     */
    private void interp_Wo(Model interp, Model prev, Model next) {
        interp_Wo2(interp, prev, next, 0.5F);
    }

    /**
     * Weighted interpolation of two Wo samples.
     *
     * @param model
     * @param prev
     * @param next
     * @param weight
     */
    private void interp_Wo2(Model model, Model prev, Model next, float weight) {
        /*
         * trap corner case where voicing est is probably wrong
         */

        if (model.getVoiced() && !prev.getVoiced() && !next.getVoiced()) {
            model.setVoiced(false);
        }

        /*
         * Wo depends on voicing of this model and adjacent models
         */
        if (model.getVoiced()) {
            if (prev.getVoiced() && next.getVoiced()) {
                model.setWo((1.0F - weight) * prev.getWo() + weight * next.getWo());
            } else if (!prev.getVoiced() && next.getVoiced()) {
                model.setWo(next.getWo());
            } else if (prev.getVoiced() && !next.getVoiced()) {
                model.setWo(prev.getWo());
            }
        } else {
            model.setWo(Defines.TWO_PI / Defines.P_MAX);
        }

        int templ = (int) ((float) Math.PI / model.getWo());    // goes wild when switching modes
        model.setL(templ > 80 ? 80 : templ);
    }

    /**
     * Interpolates center 10ms sample of energy given two samples 20ms apart.
     *
     * @param prev_e
     * @param next_e
     * @return
     */
    private float interp_energy(float prev_e, float next_e) {
        return (float) Math.pow(10.0, (Math.log10(prev_e) + Math.log10(next_e)) / 2.0);
    }

    /**
     * Interpolates center 10ms sample of energy given two samples 20ms apart.
     *
     * @param prev_e
     * @param next_e
     * @param weight
     * @return
     */
    private float interp_energy2(float prev_e, float next_e, float weight) {
        return (float) Math.pow(10.0, (1.0 - weight) * Math.log10(prev_e) + weight * Math.log10(next_e));
    }

    /**
     * Weighted interpolation of LSPs.
     *
     * @param interp
     * @param prev
     * @param next
     * @param weight
     */
    private void interpolate_lsp(float[] interp, float[] prev, float[] next, float weight) {
        for (int i = 0; i < Defines.LPC_ORD; i++) {
            interp[i] = (1.0F - weight) * prev[i] + weight * next[i];
        }
    }
}
