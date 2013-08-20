package modem;

import org.apache.commons.lang3.ArrayUtils;

import java.nio.ByteBuffer;
import java.nio.FloatBuffer;

public class Modulator {

    public static final int SAMPLE_RATE = 44100;
    public static final int SAMPLES_PER_BIT = 16;
    public static final int SAMPLES_PER_BYTE = SAMPLES_PER_BIT * 8;
    public static final int BYTES_PER_SAMPLE = 2;

    public static final byte[] SYNC_READY_BYTES;

    public static final float[] SYNC_SEQUENCE = {+1, +1, +1, -1, -1, +1, -1};
    public static final int SYNC_SAMPLES_PER_SYMBOL = 3;

    static {
        float[] tempSync = new float[SYNC_SEQUENCE.length * SYNC_SAMPLES_PER_SYMBOL];
        for (int i = 0; i < SYNC_SEQUENCE.length; i++) {
            for (int j = 0; j < SYNC_SAMPLES_PER_SYMBOL; j++) {
                tempSync[i * SYNC_SAMPLES_PER_SYMBOL + j] =  SYNC_SEQUENCE[i];
            }
        }

        SYNC_READY_BYTES = getSampleBytes(tempSync);
    }


    private static final float[] FREQUENCIES = {689, 2756};   //0 - 689 MHz, 1 -2756 MHz
    private float angle = 0;

    public byte[] getSyncBytes(){


        float samples[] = new float[64*4];
        float increment;
        for(int i = 0; i<4;i++){
            if(i==0){
                increment = (float) (2 * Math.PI) * 1378 / 44100;
            }else if(i==1){
                increment = (float) (2 * Math.PI) * 689 / 44100;
            }else if(i==2){
                increment = (float) (2 * Math.PI) * 2756 / 44100;
            }else{
                increment = (float) (2 * Math.PI) * 2067 / 44100;
            }

            for(int j=0;j<64;j++){

                samples[i*64+j] = (float) Math.sin(angle);
                angle += increment;
            }

        }


        return getSampleBytes(samples);
    }


    public byte[] modulate(byte[] data) {
        ByteBuffer samplesByteBuffer = ByteBuffer.allocate((data.length * SAMPLES_PER_BYTE * BYTES_PER_SAMPLE) + SYNC_READY_BYTES.length);
        samplesByteBuffer.put(SYNC_READY_BYTES);


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


    public static byte[] getSampleBytes(float[] samples) {
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
