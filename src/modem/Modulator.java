package modem;

import utils.Utils;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;

public class Modulator {

    public static final int SAMPLE_RATE = 44100;
    public static final int SAMPLES_PER_DATA_BIT = 16;
    public static final int SAMPLES_PER_DATA_BYTE = SAMPLES_PER_DATA_BIT * 8;
    public static final int BYTES_PER_AUDIO_SAMPLE = 2;
    public static final float MODULATION_FREQUENCY = 2756.25f;
    public static final float[] SYNC_SAMPLES;
    public static final byte[] SYNC_READY_BYTES; // Bytes of sync samples
    public static final float[] ZERO_SAMPLES;
    public static final float[] ONE_SAMPLES;

    static {
        float[] tempSync = new float[SAMPLES_PER_DATA_BIT * 2];

        float[] upHalfWave = getWave(8, 0);
        float[] downHalfWave = getWave(8, Math.PI);

        System.arraycopy(upHalfWave, 0, tempSync, 0, upHalfWave.length);
        System.arraycopy(upHalfWave, 0, tempSync, 8, upHalfWave.length);
        System.arraycopy(upHalfWave, 0, tempSync, 16, upHalfWave.length);
        System.arraycopy(downHalfWave, 0, tempSync, 24, downHalfWave.length);
        SYNC_SAMPLES = tempSync;
        SYNC_READY_BYTES = Utils.floatsToBytes(tempSync);

        ZERO_SAMPLES = getWave(SAMPLES_PER_DATA_BIT, 0);
        ONE_SAMPLES = getWave(SAMPLES_PER_DATA_BIT, Math.PI);
    }

    public byte[] modulate(byte[] data) {
        ByteBuffer samplesByteBuffer = ByteBuffer.allocate(
                (data.length * SAMPLES_PER_DATA_BYTE * BYTES_PER_AUDIO_SAMPLE) + SYNC_READY_BYTES.length);
        samplesByteBuffer.put(SYNC_READY_BYTES);

        for (byte b : data) {
            float[] floatSamples = this.getByteSamples(b);
            byte[] byteSamples = Utils.floatsToBytes(floatSamples);
            samplesByteBuffer.put(byteSamples);
        }

        return samplesByteBuffer.array();
    }

    private float[] getByteSamples(byte b) {
        FloatBuffer floatBuffer = FloatBuffer.allocate(SAMPLES_PER_DATA_BYTE);
        byte[] bits = Utils.byteToBits(b);
        for (byte bit : bits) {
            floatBuffer.put(getBitSamples(bit));
        }
        float[] buffer = floatBuffer.array();
        return buffer;
    }

    private float[] getBitSamples(byte bit) {
        if (bit == 0) {
            return ZERO_SAMPLES;
        } else {
            return ONE_SAMPLES;
        }
    }

    public static float[] getWave(int samples, double phase) {
        float[] wave = new float[samples];
        float increment = (float) (2 * Math.PI) * MODULATION_FREQUENCY / SAMPLE_RATE;
        for (int i = 0; i < wave.length; i++) {
            wave[i] = (float) (Math.sin(phase));
            phase += increment;
        }
        return wave;
    }

}
