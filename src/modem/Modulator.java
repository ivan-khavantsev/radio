package modem;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;

public class Modulator {

    public static final int SAMPLE_RATE = 44100;
    public static final int SAMPLES_PER_BIT = 16;
    public static final int SAMPLES_PER_BYTE = SAMPLES_PER_BIT * 8;
    public static final int BYTES_PER_SAMPLE = 2;
    public static final int SYNC_SAMPLES_PER_PACKET = 10;
    public static final int SYNC_BYTES_PER_PACKET = SYNC_SAMPLES_PER_PACKET * BYTES_PER_SAMPLE;
    public static final byte[] SYNC_BYTES;
    public static final float[] SYNC_PACKET = {+1, +1, +1, -1, -1, +1, -1};
    public static final int SYNC_SAMPLES_PER_SYMBOL = 3;
    public static final float[] SYNC_SAMPLES;

    static {
        float[] tempSync = new float[SYNC_PACKET.length * SYNC_SAMPLES_PER_SYMBOL];
        for (int i = 0; i < SYNC_PACKET.length; i++) {
            for (int j = 0; j < SYNC_SAMPLES_PER_SYMBOL; j++) {
                tempSync[i * SYNC_SAMPLES_PER_SYMBOL + j] =  SYNC_PACKET[i];
            }
        }
        SYNC_SAMPLES = tempSync;
        SYNC_BYTES = getSampleBytes(SYNC_SAMPLES);
    }


    private static final float[] FREQUENCIES = {689, 2756};   //0 - 689 MHz, 1 -2756 MHz
    private float angle = 0;

    public byte[] modulate(byte[] data) {
        ByteBuffer samplesByteBuffer = ByteBuffer.allocate((data.length * SAMPLES_PER_BYTE * BYTES_PER_SAMPLE) + SYNC_BYTES.length);
        samplesByteBuffer.put(SYNC_BYTES);


        for (byte b : data) {
            float[] floatSamples = this.getByteFloats(b);
            byte[] byteSamples = getSampleBytes(floatSamples);
            samplesByteBuffer.put(byteSamples);
        }


        angle = 0;
        return samplesByteBuffer.array();
    }

    private int[] getBinary(byte b) {
        int[] bits = new int[8];
        bits[0] = (b & 0x80) >> 7;
        bits[1] = (b & 0x40) >> 6;
        bits[2] = (b & 0x20) >> 5;
        bits[3] = (b & 0x10) >> 4;
        bits[4] = (b & 0x08) >> 3;
        bits[5] = (b & 0x04) >> 2;
        bits[6] = (b & 0x02) >> 1;
        bits[7] = (b & 0x01);
        return bits;
    }

    private float[] getByteFloats(byte b) {
        FloatBuffer floatBuffer = FloatBuffer.allocate(SAMPLES_PER_BYTE);
        int[] bits = this.getBinary(b);
        for (int bit : bits) {
            floatBuffer.put(freqBit(bit));
        }
        float[] buffer = floatBuffer.array();
        return buffer;
    }


    private static byte[] getSampleBytes(float[] samples) {
        byte[] buffer = new byte[samples.length * 2];
        for (int i = 0, j = 0; i < samples.length; i++, j += 2) {
            short value = (short) (samples[i] * Short.MAX_VALUE);
            buffer[j] = (byte) (value | 0xff);
            buffer[j + 1] = (byte) (value >> 8);
        }
        return buffer;
    }


    private float[] freqBit(int bit) {
        float frequency = FREQUENCIES[bit];
        float increment = (float) (2 * Math.PI) * frequency / SAMPLE_RATE;
        float samples[] = new float[SAMPLES_PER_BIT];
        for (int i = 0; i < samples.length; i++) {
            samples[i] = (float) Math.sin(angle);
            angle += increment;
        }
        return samples;
    }

}
