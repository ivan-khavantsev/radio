package modem;

import utils.Utils;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;

public class Modulator {

    public static final int SAMPLE_RATE = 44100;
    public static final int SAMPLES_PER_DATA_BIT = 16;
    public static final int SAMPLES_PER_DATA_BYTE = SAMPLES_PER_DATA_BIT * 8;
    public static final int BYTES_PER_AUDIO_SAMPLE = 2;

    public static final byte[] SYNC_READY_BYTES;
    public static final float[] SYNC_SEQUENCE = {+1, +1, +1, -1, -1, +1, -1};
    public static final int SYNC_SAMPLES_PER_SYMBOL = 4;

    private static final float[] FREQUENCIES = {689, 2756};   //0 - 689 MHz, 1 -2756 MHz
    private float angle = 0;

    static {
        float[] tempSync = new float[SYNC_SEQUENCE.length * SYNC_SAMPLES_PER_SYMBOL];
        for (int i = 0; i < SYNC_SEQUENCE.length; i++) {
            for (int j = 0; j < SYNC_SAMPLES_PER_SYMBOL; j++) {
                tempSync[i * SYNC_SAMPLES_PER_SYMBOL + j] = SYNC_SEQUENCE[i];
            }
        }
        SYNC_READY_BYTES = Utils.floatsToBytes(tempSync);
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

        angle = 0;
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
        float frequency = FREQUENCIES[bit];
        float increment = (float) (2 * Math.PI) * frequency / SAMPLE_RATE;
        float samples[] = new float[SAMPLES_PER_DATA_BIT];
        for (int i = 0; i < samples.length; i++) {
            samples[i] = (float) Math.sin(angle);
            angle += increment;
        }
        return samples;
    }

}
